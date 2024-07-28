classdef circuitMC < circuit
    % CIRCUITMC Container and controller class for circuits which include 
    % subprocesses modeled by Markov Chains.
    
    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------

    properties (SetAccess = protected)
        % internal MC variables
        Si;                  % MC states
        Ti;                  % MC mean event time
        stateChangeMatrixi;     %state change matrix for the global state vector
        stateChangeMatrixiDiff; %state change matrix for the global state vector to be used for diffusion part of SDE
        rateVectori;          % rate vectors of MCs
        dRateVectori;         % derivatives of rate vectors of MCs wrt. electricals
        dRateVectorStatei;    % derivatives of rate vectors of MCs wrt. MC states
        dQsti;               % Jacobian of charges with respect to MC states
        dIsti;               % Jacobian of currents with respect to MC states
        dRateVectorinc;       % derivatives of rate vectors of MCs wrt. electricals - no currents
        dQstinc;             % Jacobian of charges with respect to MC states - no-currents
        dIstinc;             % Jacobian of currents with respect to MC states - no-currents
        ttri = -Inf;         % Time of last time constant update
        Xtri;                % state vector - voltages - currents - MCs
        Vtri;                % no current state vector - voltages - MCs
    end
    
    methods
        
        function newCircuitMC = circuitMC(varargin)
            newCircuitMC = newCircuitMC@circuit(varargin{:});
        end
        
        function newComp = addComponent(thisCircuitTr, circComp, varargin)
            newComp = addComponent@circuit(thisCircuitTr, circComp, varargin{:});
            nMC = circComp.numMCs;
            nCircComp = thisCircuitTr.numCircCompMCs;
            if nMC > 0
                thisCircuitTr.MCs = [thisCircuitTr.MCs, circComp.MCs];
                thisCircuitTr.numMCs = thisCircuitTr.numMCs + nMC;
                
                thisCircuitTr.circCompMCs{nCircComp+1} = circComp;
                thisCircuitTr.numCircCompMCs = thisCircuitTr.numCircCompMCs + 1;
                
                thisCircuitTr.numStates = thisCircuitTr.numStates + circComp.numStates;
                thisCircuitTr.numTransitions = thisCircuitTr.numTransitions + circComp.numTransitions;
            end
        end
        
        function seal(thisCircuitTr)
            seal@circuit(thisCircuitTr);
            nMC = thisCircuitTr.numMCs;
            nV = thisCircuitTr.numVars;
            nVnc = thisCircuitTr.numVarsnc;
            nVV = thisCircuitTr.numVoltageVars;
            
            for i=1:nMC
                MC = thisCircuitTr.MCs{i};
                parent = MC.parentComponent;
                voltEqNums = parent.voltEqNums;
                currEqNums = parent.currEqNums;
                MC_EqNums = parent.MCEqNums;
                
                parent = parent.parentCircuit;
                while isequal(parent, thisCircuitTr) == 0
                    voltEqNums = parent.voltEqNums(voltEqNums);
                    currEqNums = parent.currEqNums(currEqNums);
                    MC_EqNums = parent.MCEqNums(MC_EqNums);
                    parent = parent.parentCircuit;
                end
                
                MC.eqNums = [voltEqNums, nVV + currEqNums];
                MC.neighborNums = MC_EqNums;
            end
            
            thisCircuitTr.Si = NaN(thisCircuitTr.numStates,1);
            thisCircuitTr.Xtri = NaN(nV,1);
            thisCircuitTr.Vtri = NaN(nVnc,1);
        end
        
        
        function sQ = Q(thisCircuitTr, x, t, s)
            % Q Get the charge vector Q.
            if nargin == 4
                thisCircuitTr.setMCStates(s);
            end
            sQ = Q@circuit(thisCircuitTr,x,t);
        end
        
        
        function sI = I(thisCircuitTr, x, t, s)
            % I Get the current vector I.
            if nargin == 4
                thisCircuitTr.setMCStates(s);
            end
            sI = I@circuit(thisCircuitTr,x,t);
        end
        
        
        function sJ = J(thisCircuitTr, x, t, s)
            % J Get the source vector J.
            if nargin == 4
                thisCircuitTr.setMCStates(s);
            end
            sJ = J@circuit(thisCircuitTr,x,t);
        end
        
        
        function sdQ = dQ(thisCircuitTr, x, t, s)
            % DQ Get the Jacobian of Q.
            if nargin == 4
                thisCircuitTr.setMCStates(s);
            end
            sdQ = dQ@circuit(thisCircuitTr,x,t);
        end
        
        
        function sdI = dI(thisCircuitTr, x, t, s)
            % DI Get the Jacobian of I.
            if nargin == 4
                thisCircuitTr.setMCStates(s);
            end
            sdI = dI@circuit(thisCircuitTr,x,t);
        end
        
        
        function [sQ, sI, sJ, sdQ, sdI] = stamp(thisCircuitTr, x, t, s)
            if nargin == 4
                thisCircuitTr.setMCStates(s);
            end
            [sQ, sI, sJ, sdQ, sdI] = stamp@circuit(thisCircuitTr,x,t);
        end
        
        
        function T = T(thisCircuitTr, x, t, s)
            thisCircuit.statusCheckMC(x,t,s);
            T = thisCircuitTr.Ti;
        end
        
        function rateVector = rateVector(thisCircuitTr, x, t, s)
            thisCircuitTr.statusCheckMC(x,t,s);
            rateVector = thisCircuitTr.rateVectori;
        end
        
        function dRateVector = dRateVector(thisCircuitTr, x, t, s)
            thisCircuitTr.statusCheckMC(x,t,s);
            dRateVector = thisCircuitTr.dRateVectori;
        end
        
        function dRateVectorState = dRateVectorState(thisCircuitTr, x, t, s)
            thisCircuitTr.statusCheckMC(x,t,s);
            dRateVectorState = thisCircuitTr.dRateVectorStatei;
        end
        
        
        function dQst = dQst(thisCircuitTr, x, t, s)
            % DQst Get the Jacobian of Q with respect to MC states.
            thisCircuitTr.statusCheckMC(x,t,s);
            dQst = thisCircuitTr.dQsti;
        end
        
        
        function dIst = dIst(thisCircuitTr, x, t, s)
            % DIst Get the Jacobian of I with respect to MC states.
            thisCircuitTr.statusCheckMC(x,t,s);
            dIst = thisCircuitTr.dIsti;
        end
        
        
        function [rateVector, dRateVector, dRateVectorState, dQst, dIst, stateChangeMatrix] = ttc(thisCircuitTr, x, t, s)
            
            thisCircuitTr.statusCheckMC(x,t,s);
            
            rateVector = thisCircuitTr.rateVectori;
            
            dRateVector = thisCircuitTr.dRateVectori;
            
            dRateVectorState = thisCircuitTr.dRateVectorStatei;
            
            dQst = thisCircuitTr.dQsti;
            dIst = thisCircuitTr.dIsti;
            
            stateChangeMatrix = thisCircuitTr.stateChangeMatrixi;
        end
        
        
        
        %------------------no current case--------------------------------
        
        function sQnc = Qnc(thisCircuitTr, v, t, s)
            % Q Get the charge vector Q.
            if nargin == 4
                thisCircuitTr.setMCStates(s);
            end
            sQnc = Qnc@circuit(thisCircuitTr,v,t);
        end
        
        function sInc = Inc(thisCircuitTr, v, t, s)
            % I Get the current vector I.
            if nargin == 4
                thisCircuitTr.setMCStates(s);
            end
            sInc = Inc@circuit(thisCircuitTr,v,t);
        end
        
        function sJnc = Jnc(thisCircuitTr, v, t, s)
            % J Get the source vector J.
            if nargin == 4
                thisCircuitTr.setMCStates(s);
            end
            sJnc = Jnc@circuit(thisCircuitTr,v,t);
        end
        
        function sdQnc = dQnc(thisCircuitTr, v, t, s)
            % DQ Get the Jacobian of Q.
            if nargin == 4
                thisCircuitTr.setMCStates(s);
            end
            sdQnc = dQnc@circuit(thisCircuitTr,v,t);
        end
        
        function sdInc = dInc(thisCircuitTr, v, t, s)
            % DI Get the Jacobian of I.
            if nargin == 4
                thisCircuitTr.setMCStates(s);
            end
            sdInc = dInc@circuit(thisCircuitTr,v,t);
        end
        
        function [sQ, sI, sJ, sdQ, sdI] = stampnc(thisCircuitTr, v, t, s)
            if nargin == 4
                thisCircuitTr.setMCStates(s);
            end
            [sQ, sI, sJ, sdQ, sdI] = stampnc@circuit(thisCircuitTr,v,t);
        end
        
        % these function will give the same results as the ones above.
        
        function Tnc = Tnc(thisCircuitTr, v, t, s)
            thisCircuitTr.statusCheckNoCurrMC(v,t,s);
            Tnc = thisCircuitTr.Ti;
        end
        
        function rateVectornc = rateVectornc(thisCircuitTr, v, t, s)
            thisCircuitTr.statusCheckNoCurrMC(v,t,s);
            rateVectornc = thisCircuitTr.rateVectori;
        end
        
        function dRateVectorStatenc = dRateVectorStatenc(thisCircuitTr, v, t, s)
            thisCircuitTr.statusCheckNoCurrMC(v,t,s);
            dRateVectorStatenc = thisCircuitTr.dRateVectorStatei;
        end
        
        function dRateVectornc = dRateVectornc(thisCircuitTr, v, t, s)
            thisCircuitTr.statusCheckNoCurrMC(v,t,s);
            dRateVectornc = thisCircuitTr.dRateVectorinc;
        end
        
        
        function dQstnc = dQstnc(thisCircuitTr, v, t, s)
            % DQstnc Get the Jacobian of Q with respect to MC states.
            % - no current variables
            thisCircuitTr.statusCheckNoCurrMC(v,t,s);
            dQstnc = thisCircuitTr.dQstinc;
        end
        
        
        function dIstnc = dIstnc(thisCircuitTr, v, t, s)
            % DIstnc Get the Jacobian of I with respect to MC states.
            % - no current variables.
            thisCircuitTr.statusCheckNoCurrMC(v,t,s);
            dIstnc = thisCircuitTr.dIstinc;
        end
        
        
        function [rateVector, dRateVectornc, dRateVectorState, dQstnc, dIstnc, stateChangeMatrix] = ...
                ttcnc(thisCircuitTr, v, t, s)
            thisCircuitTr.statusCheckNoCurrMC(v,t,s);
            
            rateVector = thisCircuitTr.rateVectori;

            dRateVectornc = thisCircuitTr.dRateVectorinc;
            
            dRateVectorState = thisCircuitTr.dRateVectorStatei;
            
            dQstnc = thisCircuitTr.dQstinc;
            dIstnc = thisCircuitTr.dIstinc;
            
            stateChangeMatrix = thisCircuitTr.stateChangeMatrixi;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %                    MC handling                                %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function setMCStates(thisCircuitTr, MCStates)
            if length(MCStates) ~= thisCircuitTr.numStates
                error('State vector does not have the correct size!');
            end
            k=0;
            for i=1:thisCircuitTr.numMCs
                m=thisCircuitTr.MCs{i}.numStates;
                thisCircuitTr.MCs{i}.setState(MCStates(k+1:k+m,1));
                k=k+m;
            end
        end
        
        
        % These functions all assume that the internal state of the circuit
        % is up-to-date. We don't want to use ttc functions here that
        % calculate time constants based on a circuit state that is
        % consistent with the MC states because at a MC firing time the
        % voltages stay constant although the MC state changes.
        
        function Lambda = MCEventRates(thisCircuitTr, x, t, s)
            [rateVector] = thisCircuitTr.ttc(x,t,s);
            Lambda = rateVector;
        end
        
        function Lambda = MCEventRatesNoCurr(thisCircuitTr, v, t, s)
            [rateVector] = thisCircuitTr.ttcnc(v,t,s);
            if isempty(rateVector)
                %[m,n]=size(rateVector);
                m = thisCircuitTr.numTransitions;
                Lambda = Inf*ones(m,1);
            else
                Lambda = rateVector;
            end
        end
        
        function [rtot, rs] = getTotalMCEventRate(thisCircuitTr)
            % rtot ... total rate
            % rs ... individual rates
            % assume time constants are up-to-date
            S = thisCircuit.Si;
            nTr = thisCircuitTr.numMCs;
            
            for i=1:1:nTr
                S_floor=floor(S);
                
                if any(S_floor ~= S)
                    error('MC states must be integer!');
                end
            end
            
            rs=thisCircuit.rateVectori;
            
            rtot = sum(rs);
        end
        
        function flipMCState(thisCircuitTr, indMC, indEvent)
            thisCircuitTr.MCs{indMC}.flipState(indEvent);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %                    MC handling  OLD                           %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function ts = getMCStates(thisCircuitTr)
            ts = zeros(thisCircuitTr.numTotalStates,1);
            
            k=0;
            for i = 1:thisCircuitTr.numMCs
                m = thisCircuitTr.MCs{i}.numStates;
                ts(k+1:k+m,1) = thisCircuitTr.MCs{i}.stateVector;
                k=k+m;
            end
        end
    end
    
    methods (Access = protected)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %               MC related equation generation                  %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function statusCheckMC(thisCircuitTr, x, t, s)
            % STATUSCHECK check if the internal state matches the new
            % state (x,t).
            
            if length(s) ~= thisCircuitTr.numStates
                error('The MC state vector does not have the correct size!');
            end
            
            k=0;
            for i=1:thisCircuitTr.numMCs
                m = thisCircuitTr.MCs{i}.numStates;
                if any(thisCircuitTr.MCs{i}.stateVector ~= s(k+1:k+m,1))
                    % this will set forceEquation flag and statusCheck will
                    % recompute equations.
                    thisCircuitTr.setMCStates(s);
                    break; %there is no need to further check for difference of the states
                end
            end
            
            status = thisCircuitTr.statusCheck(x,t);
            if status || (thisCircuitTr.ttri ~= t) || any(thisCircuitTr.Xtri ~= x)
                thisCircuitTr.ttcEval(s);
            end
        end
        
        function statusCheckNoCurrMC(thisCircuitTr, v, t, s)
            % number of independent variables in the equations are set as
            % (number of independent voltage variables) - (number of
            % voltage sources) because one node gets eliminated per voltage
            % source and treated as a dependent to the remaining node.
            % status is true if we recompute the equations
            
            if length(s) ~= thisCircuitTr.numStates
                error('The MC state vector does not have the correct size!');
            end
            
            k=0;
            for i=1:thisCircuitTr.numMCs
                m = thisCircuitTr.MCs{i}.numStates;
                if any(thisCircuitTr.MCs{i}.stateVector ~= s(k+1:k+m,1))
                    % this will set forceEquation flag and statusCheck will
                    % recompute equations.
                    thisCircuitTr.setMCStates(s);
                    break; %there is no need to further check for difference of the states
                end
            end
            
            status = thisCircuitTr.statusCheckNoCurr(v,t);
            if status || (thisCircuitTr.ttri ~= t) || any(thisCircuitTr.Vtri ~= v)
                thisCircuitTr.ttcEval(s);
            end
        end
        
        function ttcEval(thisCircuitTr, s)
            nEq = thisCircuitTr.numVoltageVars + thisCircuitTr.numCurrentVars;
            eqNumsnc = thisCircuitTr.eqNumsnc;
            gndNum = thisCircuitTr.groundNodeNumber;
            
            nCircCompMC = thisCircuitTr.numCircCompMCs;
            nStates = thisCircuitTr.numStates;
            nTransitions = thisCircuitTr.numTransitions;
            
            rateVector = zeros(nTransitions,1);
            
            stateChangeMatrix = [];
            stateChangeMatrixDiff = [];
            
            dRateVector = zeros(nTransitions, nEq);
            
            dRateVectorState = zeros(nTransitions, nStates);
            
            dQst = zeros(nEq, nStates);
            dIst = zeros(nEq, nStates);
            
            i=1;
            m=0;
            n=0;
            k=m;
            l=n;
            nTr = thisCircuitTr.numMCs;
            parentComponentCurrent = thisCircuitTr.MCs{1}.parentComponent;
            statesCircComp = parentComponentCurrent.numStates;
            transitionsCircComp = parentComponentCurrent.numTransitions;
            
            for j=1:nTr
                
                tr = thisCircuitTr.MCs{j};
                
                if isequal(parentComponentCurrent, tr.parentComponent) == false
                    parentComponentPrev = parentComponentCurrent;
                    parentComponentCurrent = tr.parentComponent;
                    statesCircComp = parentComponentCurrent.numStates;
                    transitionsCircComp = parentComponentCurrent.numTransitions;
                    m = m + parentComponentPrev.numTransitions;
                    n = n + parentComponentPrev.numStates;
                    k = m;
                    l = n;
                    i = i+1;
                end
                
                trEqNums = tr.neighborNums;
                eqNums = tr.eqNums;
                statesMC = tr.numStates;
                transitionsMC = tr.numTransitions;
                
                tr.updateRates();
                [sRateVector, sdRateVector, sdRateVectorState, sdQst, sdIst] = tr.getRates();
                
                rateVector(k+1:k+transitionsMC,1) = sRateVector;
                stateChangeMatrix(l+1:l+statesMC,k+1:k+transitionsMC) = tr.stateChangeMatrix;
                
                if strcmp(tr.parentComponent.type, 'HH') == true
                    stateChangeMatrixDiff(l+1:l+statesMC,k+1:k+transitionsMC) = zeros(statesMC,transitionsMC);
                else
                    stateChangeMatrixDiff(l+1:l+statesMC,k+1:k+transitionsMC) = tr.stateChangeMatrix;
                end
                    
                dRateVector(k+1:k+transitionsMC,eqNums) = sdRateVector;
                dRateVectorState(k+1:k+transitionsMC,n+1:n+statesCircComp) = sdRateVectorState;
                dQst(eqNums,n+1:n+statesCircComp) = sdQst;
                dIst(eqNums,n+1:n+statesCircComp) = sdIst;
                
                if isa(tr,'inReceptorMC') == true || isa(tr,'exReceptorMC') == true || isa(tr,'exSlowReceptorMC') == true
                    sourceNeuron = tr.sourceNeuron;
                    sourceEqNums = sourceNeuron.voltEqNums;
                    [sdRateVectorCross]=tr.getCrossRates();
                    dRateVector(k+1:k+transitionsMC,sourceEqNums) = dRateVector(k+1:k+transitionsMC,sourceEqNums) + sdRateVectorCross;
                end
                k=k+transitionsMC;
                l=l+statesMC;
            end
            
            thisCircuitTr.Si = s;
            thisCircuitTr.stateChangeMatrixi = stateChangeMatrix;
            thisCircuitTr.stateChangeMatrixiDiff = stateChangeMatrixDiff;
            
            thisCircuitTr.Vtri = thisCircuitTr.Vi;
            thisCircuitTr.Xtri = thisCircuitTr.Xi;
            thisCircuitTr.ttri = thisCircuitTr.ti;
            
            thisCircuitTr.rateVectori = rateVector;
            thisCircuitTr.dRateVectorStatei = dRateVectorState;
            
            if gndNum ~= 0
                idx = [1:gndNum-1, gndNum+1:nEq];
                
                thisCircuitTr.dRateVectori=dRateVector(:,idx);
                thisCircuitTr.dQsti=dQst(idx,:);
                thisCircuitTr.dIsti=dIst(idx,:);
            else
                thisCircuitTr.dRateVectori = dRateVector;
                thisCircuitTr.dQsti=dQst;
                thisCircuitTr.dIsti=dIst;
            end
            
            % current free case
            thisCircuitTr.dRateVectorinc = dRateVector(:,eqNumsnc);
            thisCircuitTr.dQstinc = dQst(eqNumsnc, :);
            thisCircuitTr.dIstinc = dIst(eqNumsnc, :);
        end
         
    end
end

