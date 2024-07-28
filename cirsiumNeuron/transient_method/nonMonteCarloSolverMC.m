classdef nonMonteCarloSolverMC < handle
    % NONMONTECARLOSOLVERMC time-domain, semi-analytical, non Monte Carlo
    % solver for a neuronal circuit. This solver is used for transient
    % stochastic characterization of a neuronal circuit. Correlation matrix
    % for circuit variables is computed.
    
    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------
    
    properties (Constant)
        TEROH = 0.1; % MC event rate overhead
        
        keyNames = {'solverType'
            'linSolverType'
            'stateSpaceType'
            'currentVariables'
            't0'
            'tend'
            'timeStep'
            'breakPoints'
            'reltol'
            'abstol_v'
            'abstol_c'
            'abstol_q'
            'maxOrder'
            'maxStep'
            'chargeScaleFactor'
            'y0'
            'yp0'
            's0'
            'sp0'
            'seq0'
            'showSolution'
            'lmax'};
        validSolverTypes = {'TRPZ-BE', 'TRPZ-RK'};
        validLinSolverTypes = {'Dense', 'GMRES'};
        validStateSpaceTypes = {'reduced', 'standard'}
        validOnOff = {'on', 'off'};
    end
    
    
    properties (SetAccess = private)
        circuit = circuit.empty;
        solverType = 'TRPZ-RK';
        linSolverType = 'GMRES';
        stateSpaceType = 'standard';
        currentVariables = 'off';
        showSolution = 'on';
        resFunIC;
        jacFunIC;
        precSolveFunIC;
        resFun;
        jacFun;
        precSolveFun;
        QFun;
        IFun;
        JFun;
        dQFun;
        dIFun;
        dQstFun;
        dIstFun;
        stampFun;
        ttcFun;
        t0 = 0;
        breakPoints;
        reltol = 1e-3;
        abstol_v = 1e-6;
        abstol_c = 1e-9;
        abstol_q = 1e-21;
        restol = 1e-6;
        tolVecIC;
        maxOrder = 5;
        maxStep = inf;
        timeStep = 1e-12;
        chargeScaleFactor = 1e9;
        y0;
        yp0;
        y0I;
        yp0I;
        s0;
        sp0;
        seq0;
        seqp0;
        ynow;
        ypnow;
        Qprev;
        tprev;
        tnow;
        tref;
        tend;
        tolVec;
        IDAOptions;
        Y;
        T;
        Tref;
        SC;
        SV;
        
        PROPENSITYVECTOR;
        RATEVECTOR;
        
        stats;
        numVars;
        numVoltageVars;
        numNewtonIters;
        condNumbers;
        nextStopTime;
        nextBreakPoint;
        
        numMCs;
        
        numStates;
        numTransitions;
        
        snow;
        lmax;
        lnow;
        
        taucnow;
        tauenow;
        propensityVectorNow;
        rateVectorNow;
        
        loverhead = nonMonteCarloSolverMC.TEROH;
        S;
        Seq;
        K;
        Kcur;
        Kprev;
        
        Kref;
        Krefcur;
        Krefprev;
        Kind;
        eventTimes;
        nextEventTime;
        
        firedMCs;
        
        numBreaks = 0;
        numIters = 0;
        initialDerivative;
        calcCurrents = false;
        lastState;
    end
    
    methods
        
        function obj = nonMonteCarloSolverMC(circ, varargin)
            if nargin > 0
                if nargin > 1
                    obj.setOptions(varargin{:});
                end
                obj.setCircuit(circ);
            end
        end
        
        function setOptions(thisSolver, varargin)
            if rem(nargin-1, 2) ~= 0
                error('Options structure is missing a key or a value!');
            end
            for i=1:2:nargin-1
                key = varargin{i};
                if ~ischar(key)
                    error('Options must consist of key-value pairs');
                end
                val = varargin{i + 1};
                try
                    thisSolver.(key) = val;
                catch e
                    disp(e);
                    error('%s is an unrecognized key!', key);
                end
            end
            
            % check if the options are valid
            
            if all(strcmp(thisSolver.solverType, nonMonteCarloSolverMC.validSolverTypes) == 0)
                msgStr = sprintf('Unrecognized solver type: %s', thisSolver.solverType);
                msgStr = sprintf('%s\nValid choices are:', msgStr);
                for i=1:numel(nonMonteCarloSolverMC.validSolverTypes)
                    msgStr = sprintf('%s\n%d. %s', msgStr, i, nonMonteCarloSolverMC.validSolverTypes{i});
                end
                error(msgStr);
            end
            
            if all(strcmp(thisSolver.linSolverType, nonMonteCarloSolverMC.validLinSolverTypes) == 0)
                msgStr = sprintf('Unrecognized linear solver type: %s', thisSolver.linSolverType);
                msgStr = sprintf('%s\nValid choices are:', msgStr);
                for i=1:numel(nonMonteCarloSolverMC.validLinSolverTypes)
                    msgStr = sprintf('%s\n%d. %s', msgStr, i, nonMonteCarloSolverMC.validLinSolverTypes{i});
                end
                error(msgStr);
            end
            
            if all(strcmp(thisSolver.stateSpaceType, nonMonteCarloSolverMC.validStateSpaceTypes) == 0)
                msgStr = sprintf('Unrecognized state space type: %s', thisSolver.stateSpaceType);
                msgStr = sprintf('%s\nValid choices are:', msgStr);
                for i=1:numel(nonMonteCarloSolverMC.validStateSpaceTypes)
                    msgStr = sprintf('%s\n%d. %s', msgStr, i, nonMonteCarloSolverMC.validStateSpaceTypes{i});
                end
                error(msgStr);
            end
            
            if all(strcmp(thisSolver.currentVariables, nonMonteCarloSolverMC.validOnOff) == 0)
                msgStr = sprintf('Unrecognized current variables option: %s', thisSolver.currentVariables);
                msgStr = sprintf('%s\nValid choices are:', msgStr);
                for i=1:numel(nonMonteCarloSolverMC.validOnOff)
                    msgStr = sprintf('%s\n%d. %s', msgStr, i, nonMonteCarloSolverMC.validOnOff{i});
                end
                error(msgStr);
            end
            
            if all(strcmp(thisSolver.showSolution, nonMonteCarloSolverMC.validOnOff) == 0)
                msgStr = sprintf('Unrecognized plotting option: %s', thisSolver.showSolution);
                msgStr = sprintf('%s\nValid choices are:', msgStr);
                for i=1:numel(nonMonteCarloSolverMC.validOnOff)
                    msgStr = sprintf('%s\n%d. %s', msgStr, i, nonMonteCarloSolverMC.validOnOff{i});
                end
                error(msgStr);
            end
            
            % did we get an initial derivative vector?
            thisSolver.initialDerivative = ~isempty(thisSolver.yp0);
            
            % if so we can calculate a steady state solution etc.
            
            % set the solver functions for initial condition calculations
            % with IDA
            switch thisSolver.linSolverType
                case 'GMRES'
                    jacFunName = 'IDAJacTimesVecFun';
                case 'Dense'
                    jacFunName = 'IDADJacFun';
            end
            
            
            if strcmp(thisSolver.currentVariables, 'off')
                thisSolver.calcCurrents = true;
            end
            
            resFunName = 'IDAResFun';
            precSolveFunName = 'IDAPrecSolveFun';
            if strcmp(thisSolver.solverType, 'TRPZ-BE') || strcmp(thisSolver.solverType, 'TRPZ-RK') || strcmp(thisSolver.solverType, 'IDA-BE')
                resFunName = [resFunName, 'EMS'];
                jacFunName = [jacFunName, 'EMS'];
                precSolveFunName = [precSolveFunName, 'EMS'];
            end
            
            thisSolver.resFunIC = eval(['@thisSolver.' resFunName]);
            thisSolver.jacFunIC = eval(['@thisSolver.' jacFunName]);
            thisSolver.precSolveFunIC = eval(['@thisSolver.' precSolveFunName]);
            
            % set the solver functions
            switch thisSolver.solverType
                case {'TRPZ-BE', 'TRPZ-RK'}
                    thisSolver.resFun = @thisSolver.TRPZResFun;
                    thisSolver.jacFun = @thisSolver.TRPZJacFun;
                case 'IDA-BE'
                    thisSolver.resFun = thisSolver.resFunIC;
                    thisSolver.jacFun = thisSolver.jacFunIC;
                    thisSolver.precSolveFun = thisSolver.precSolveFunIC;
            end
            
        end
        
        function setCircuit(thisSolver, circ)
            if isa(circ, 'circuit') || isa(circ, 'subcircuit')
                thisSolver.circuit = circ;
            else
                error('The solver only accepts objects of type circuit or subcircuit');
            end
            
            % are the source currents included as variables?
            switch thisSolver.currentVariables
                case 'on'
                    thisSolver.numVars = circ.numVars;
                    thisSolver.numVoltageVars = circ.numIndepVoltVars;
                case 'off'
                    thisSolver.numVars = circ.numVarsnc;
                    thisSolver.numVoltageVars = thisSolver.numVars;
            end
            
            % set the number of MCs
            thisSolver.numMCs = circ.numMCs;
            thisSolver.numStates = circ.numStates;
            thisSolver.numTransitions = circ.numTransitions;
            
            if thisSolver.initialDerivative && isempty(thisSolver.sp0)
                thisSolver.sp0 = zeros(thisSolver.numStates,1);
            end
            
            
            % check if initial conditions have the correct size
            
            if ~isempty(thisSolver.s0) && (thisSolver.numStates ~= length(thisSolver.s0))
                error('Size of s0 is not consistent with the number of MCs!');
            end
            
            if ~isempty(thisSolver.seq0) && (thisSolver.numStates ~= length(thisSolver.seq0))
                error('Size of seq0 is not consistent with the number of MCs!');
            end
            
            if isempty(thisSolver.y0)
                error(['please supply an initial condition vector or  an'...
                    'initial guess for the steady state solution']);
            elseif length(thisSolver.y0) ~= thisSolver.numVars
                error('Size of y0 is not consistent with the problem size!');
            end
            
            if (thisSolver.initialDerivative) && (length(thisSolver.yp0) ~= thisSolver.numVars)
                error('Size of yp0 is not consistent with the problem size!');
            end
            
            % construct tolerance vectors
            nV = thisSolver.numVars;
            nVV = thisSolver.numVoltageVars;
            nMCs = thisSolver.numMCs;
            nStates = thisSolver.numStates;
            Kchrg = thisSolver.chargeScaleFactor;
            
            thisSolver.tolVecIC = [ones(nVV,1).*thisSolver.abstol_v;
                ones(nV-nVV,1).*thisSolver.abstol_c;
                ones(nV,1)*thisSolver.abstol_q*thisSolver.chargeScaleFactor];
            
            switch thisSolver.solverType
                case {'TRPZ-BE', 'TRPZ-RK'}
                    thisSolver.tolVecIC = [ones(nVV,1).*thisSolver.abstol_v;
                        ones(nV-nVV,1).*thisSolver.abstol_c;
                        ones(nV,1)*thisSolver.abstol_q*Kchrg;
                        ones(nStates,1)*thisSolver.abstol_v];
                    thisSolver.tolVec = thisSolver.tolVecIC([1:nV, 2*nV+1:2*nV+nStates]);
                    thisSolver.restol = [ones(nVV,1)*thisSolver.abstol_q*Kchrg;
                        ones(nV-nVV+nStates,1)*thisSolver.abstol_v];
                case 'IDA-BE'
                    thisSolver.tolVecIC = [ones(nVV,1).*thisSolver.abstol_v;
                        ones(nV-nVV,1).*thisSolver.abstol_c;
                        ones(nV,1)*thisSolver.abstol_q*Kchrg;
                        ones(nStates,1)*thisSolver.abstol_v];
                    thisSolver.tolVec = thisSolver.tolVecIC;
            end
            
            % set native residual and jacobian functions
            switch thisSolver.currentVariables
                case 'on'
                    thisSolver.QFun = @circ.Q;
                    thisSolver.IFun = @circ.I;
                    thisSolver.JFun = @circ.J;
                    thisSolver.dQFun = @circ.dQ;
                    thisSolver.dIFun = @circ.dI;
                    thisSolver.dQstFun = @circ.dQst;
                    thisSolver.dIstFun = @circ.dIst;
                    thisSolver.stampFun = @circ.stamp;
                    if nMCs > 0
                        thisSolver.ttcFun = @circ.ttc;
                    end
                case 'off'
                    thisSolver.QFun = @circ.Qnc;
                    thisSolver.IFun = @circ.Inc;
                    thisSolver.JFun = @circ.Jnc;
                    thisSolver.dQFun = @circ.dQnc;
                    thisSolver.dIFun = @circ.dInc;
                    thisSolver.dQstFun = @circ.dQstnc;
                    thisSolver.dIstFun = @circ.dIstnc;
                    thisSolver.stampFun = @circ.stampnc;
                    if nMCs > 0
                        thisSolver.ttcFun = @circ.ttcnc;
                    end
            end
            
            
        end
        
        function setIDAOptions(thisSolver)
            nV = thisSolver.numVars;
            nStates = thisSolver.numStates;
            
            diffInd = [zeros(nV,1); ones(nV,1)];
            
            if strcmp(thisSolver.solverType, 'TRPZ-BE') || strcmp(thisSolver.solverType, 'TRPZ-RK') || strcmp(thisSolver.solverType, 'IDA-BE')
                diffInd = [diffInd; ones(nStates,1)];
            end
            
            % set the end of times as the last break point
            thisSolver.breakPoints = [thisSolver.breakPoints,...
                thisSolver.tend];
            % set the first break point as the next break point
            thisSolver.nextBreakPoint = thisSolver.breakPoints(1);
            
            % set the solver functions in the IDAOptions structure
            thisSolver.IDAOptions = IDASetOptions(...
                'LinearSolver', thisSolver.linSolverType,...
                'PrecModule', 'UserDefined',...
                'VariableTypes', diffInd,...
                'RelTol', thisSolver.reltol,...
                'AbsTol', thisSolver.tolVecIC,...
                'MaxNumSteps', 500,...
                'MaxStep', thisSolver.maxStep,...
                'MaxOrder', thisSolver.maxOrder,...
                'JacobianFn', thisSolver.jacFunIC,...
                'PrecSolveFn', thisSolver.precSolveFunIC,...
                'StopTime', thisSolver.nextBreakPoint);
        end
        
        
        function solve(thisSolver)
            % this function calls the main solver loop.
            if isempty(thisSolver.circuit)
                error(['This solver does not have a circuit to solve! ',...
                    'Please add a circuit first.']);
            elseif thisSolver.circuit.isSealed == false
                error('The circuit is not sealed! Please seal it first.');
            end
            
            thisSolver.setIDAOptions();
            thisSolver.calcIC();
            
            switch thisSolver.solverType
                case {'TRPZ-BE', 'TRPZ-RK'}
                    thisSolver.solveTRPZ;
                case 'IDA-BE'
                    thisSolver.solveIDA;
            end
            
            if strcmp(thisSolver.showSolution, 'on')
                thisSolver.displaySolution;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %           Consistent initial conditions calculation             %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function calcIC(thisSolver)
            y0 = thisSolver.y0;
            yp0 = thisSolver.yp0;
            t0 = thisSolver.t0;
            sp0 = thisSolver.sp0;
            Kchrg = thisSolver.chargeScaleFactor;
            
            % first set the initial MC states
            if thisSolver.numMCs > 0
                if isempty(thisSolver.s0)
                    thisSolver.s0 = zeros(thisSolver.numStates,1);
                    
                    k=0;
                    for i=1:1:thisSolver.numMCs
                        m = thisSolver.circuit.MCs{i}.numStates;
                        thisSolver.s0(k+1:k+m,1) = thisSolver.circuit.MCs{i}.stateVector;
                        k=k+m;
                    end
                end
                thisSolver.circuit.setMCStates(thisSolver.s0);
                s0 = thisSolver.s0;
            end
            
            y02 = [y0; Kchrg*thisSolver.QFun(y0,t0)];
            if thisSolver.initialDerivative
                % user wants the solution with specified derivatives,
                % e.g. steady state solution
                yp02 = [yp0; Kchrg*(thisSolver.dQFun(y0,t0)*yp0)];
            else
                % we are only given the algebraic initial conditions
                yp02 = [zeros(thisSolver.numVars,1);
                    Kchrg*(thisSolver.JFun(y0,t0) - thisSolver.IFun(y0,t0))];
            end
            
            if strcmp(thisSolver.solverType, 'TRPZ-BE') || strcmp(thisSolver.solverType, 'TRPZ-RK') || strcmp(thisSolver.solverType, 'IDA-BE')
                y02 = [y02; s0];
                if thisSolver.initialDerivative
                    % user wants the solution with spcified derivatives,
                    % e.g. steady state solution
                    yp02 = [yp02; sp0];
                else
                    % we are only given the algebraic initial conditions
                    yp02 = [yp02;
                        zeros(thisSolver.numStates,1)];
                end
            end
            y02 = full(y02);
            yp02 = full(yp02);
            
            IDAInit(thisSolver.resFunIC, t0, y02, yp02, thisSolver.IDAOptions);
            
            t1 = thisSolver.t0 + thisSolver.timeStep;
            
            tic
            if thisSolver.initialDerivative
                [~, thisSolver.y0I, thisSolver.yp0I] = IDACalcIC(t1, 'FindAll');
            else
                [~, thisSolver.y0I, thisSolver.yp0I] = IDACalcIC(t1, 'FindAlgebraic');
            end
            
            
            fprintf(1,'%g seconds for consistent initial conditions.\n', toc);
            
            % Free the IDA memory if we will not use it later
            if ~strcmp(thisSolver.solverType, 'IDA-BE')
                IDAFree;
            end
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %           Non Monte Carlo Time Domain TRPZ Solution               %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function solveTRPZ(thisSolver)
            
            circ = thisSolver.circuit;
            
            h = thisSolver.timeStep;
            nV = thisSolver.numVars;
            nTr = thisSolver.numMCs;
            nStates = thisSolver.numStates;
            nTransitions = thisSolver.numTransitions;
            stateChangeMatrix = circ.stateChangeMatrixi;
            stateChangeMatrixDiff = circ.stateChangeMatrixiDiff;
            
            excludedState = 1;
            excludeStateInd(1) = excludedState;
            for ind=2:1:nTr
                m = thisSolver.circuit.MCs{ind-1}.numStates;
                excludeStateInd(ind) = excludeStateInd(ind-1) + m;
            end
            
            y0I = thisSolver.y0I;
            yp0I = thisSolver.yp0I;
            Kchrg = thisSolver.chargeScaleFactor;
            
            thisSolver.T = thisSolver.t0:h:thisSolver.tend;
            thisSolver.Tref = thisSolver.tref:h:thisSolver.tend;
            numT = numel(thisSolver.T);
            numTref = numel(thisSolver.Tref);
            t = thisSolver.t0;
            tref = thisSolver.tref;
            thisSolver.tnow = thisSolver.t0;
            Kind = thisSolver.Kind;
            
            thisSolver.Y = zeros(nV,numT);
            v = y0I(1:nV);
            thisSolver.ynow = v;
            thisSolver.Y(:,1) = v;
            
            thisSolver.S = zeros(nStates, numT);
            s = y0I(2*nV+1:2*nV+nStates);
            thisSolver.snow = s;
            thisSolver.S(:,1) = s;
            
            if strcmp(thisSolver.stateSpaceType, 'standard')
                thisSolver.K = zeros(length(Kind), numT);
                thisSolver.Kref = zeros(length(Kind), numTref);
                thisSolver.Kcur = zeros(nV+nStates, nV+nStates);
                thisSolver.Krefcur = zeros(nV+nStates, nV+nStates);
                thisSolver.Kprev = zeros(nV+nStates, nV+nStates);
                thisSolver.Krefprev = zeros(nV+nStates, nV+nStates);
            elseif strcmp(thisSolver.stateSpaceType, 'reduced')
                thisSolver.K = zeros(length(Kind), numT);
                thisSolver.Kref = zeros(length(Kind), numTref);
                thisSolver.Kcur = zeros(nV+nStates-nTr, nV+nStates-nTr);
                thisSolver.Krefcur = zeros(nV+nStates-nTr, nV+nStates-nTr);
                thisSolver.Kprev = zeros(nV+nStates-nTr, nV+nStates-nTr);
                thisSolver.Krefprev = zeros(nV+nStates-nTr, nV+nStates-nTr);
            end
            
            thisSolver.PROPENSITYVECTOR = zeros(nTransitions, numT);
            
            [Qprev, Iprev, Jprev, dQprev, dIprev] = thisSolver.stampFun(v,t,s);
            
            [rateVector, dRateVector, dRateVectorState, dQst, dIst] = ...
                thisSolver.ttcFun(v,t,s);
            
            [propensityVector, dPropensityVector, dPropensityVectorStatePart1, dPropensityVectorStatePart2] = ...
                thisSolver.propensityFun(rateVector, dRateVector, dRateVectorState, s);
            dPropensityVectorState =  dPropensityVectorStatePart1 + dPropensityVectorStatePart2;
            
            %initial matrices for non Monte Carlo analysis
            Cprev = [dQprev, dQst;
                zeros(nStates,nV), eye(nStates)];
            
            if strcmp(thisSolver.stateSpaceType, 'reduced')
                stateInd=0;
                for j=1:nTr
                    m = thisSolver.circuit.MCs{j}.numStates;
                    excludeColumn = -Cprev(:,excludeStateInd(j)+nV);
                    insertMatrix = repmat(excludeColumn,1,m);
                    Cprev(:,nV+stateInd+1:nV+stateInd+m) = Cprev(:,nV+stateInd+1:nV+stateInd+m) + insertMatrix;
                    stateInd = stateInd + m;
                end
                Cprev(:,excludeStateInd+nV)=[];
                Cprev(excludeStateInd+nV,:)=[];
            end
            
            thisSolver.numNewtonIters = zeros(1,numT);
            calcCurrents = thisSolver.calcCurrents;
            
            if calcCurrents
                thisSolver.SC = zeros(circ.numVoltageSources,numT);
                thisSolver.SV = zeros(circ.numVoltageSources,numT);
                [SVt, SCt] = circ.getSourceVariables(y0I(1:nV), t, [], [],...
                    yp0I(nV+1:2*nV));
                thisSolver.SC(:,1) = SCt;
                thisSolver.SV(:,1) = SVt;
            end
            
            thisSolver.displayTime();
            
            %%first prediction for all variables
            % predictor calculations
            sprdct = s + yp0I(2*nV+1:2*nV+nStates)*h;
            vprdct = v + yp0I(1:nV)*h;
            % end predictor
            xprdct = [vprdct; sprdct];
            %xprdct = [vprdct; sprdct; seqprdct];
            
            thisSolver.PROPENSITYVECTOR(:,1) = propensityVector;
            
            for i=2:numT
                if calcCurrents == true
                    QprevFull = thisSolver.circuit.getFullQ(v(1:nV), thisSolver.tnow);
                end
                
                t = thisSolver.T(i);
                thisSolver.tnow = t;
                
                try
                    [x, numIter] = thisSolver.newtSol(thisSolver.resFun(t, Qprev, Iprev, Jprev, h, s, propensityVector),...
                        thisSolver.jacFun(t, h), xprdct);
                catch exception
                    disp('Newton solver failed');
                    error(getReport(exception));
                    break;
                end
                
                v = x(1:nV);
                s = x(nV+1:nV+nStates);
                
                %extract the charges and time constants before changing MC states
                [Qprev, Iprev, Jprev, dQprev, dIprev] = thisSolver.stampFun(v, thisSolver.tnow);
                
                [rateVector, dRateVector, dRateVectorState, dQst, dIst] = ...
                    thisSolver.ttcFun(v,t,s);
                [propensityVector, dPropensityVector, dPropensityVectorStatePart1, dPropensityVectorStatePart2] = ...
                    thisSolver.propensityFun(rateVector, dRateVector, dRateVectorState, s);
                dPropensityVectorState =  dPropensityVectorStatePart1 + dPropensityVectorStatePart2;
                
                %form the matrices for Lyapunov equation of the variance-covariance function
                Cnew = [dQprev, dQst;
                    zeros(nStates,nV), eye(nStates)];
                
                Gnew = [dIprev, dIst;
                    -stateChangeMatrix*dPropensityVector, -stateChangeMatrix*dPropensityVectorState];
                
                for j=1:1:nStates
                    B_lower(j,:) = stateChangeMatrixDiff(j,:).*sqrt(propensityVector');
                end
                Bnew = [zeros(nV, nTransitions); B_lower];
                
                if strcmp(thisSolver.stateSpaceType, 'reduced')
                    stateInd=0;
                    for j=1:nTr
                        m = thisSolver.circuit.MCs{j}.numStates;
                        excludeColumnC = -Cnew(:,excludeStateInd(j)+nV);
                        excludeColumnG = -Gnew(:,excludeStateInd(j)+nV);
                        insertMatrixC = repmat(excludeColumnC,1,m);
                        insertMatrixG = repmat(excludeColumnG,1,m);
                        Cnew(:,nV+stateInd+1:nV+stateInd+m) = Cnew(:,nV+stateInd+1:nV+stateInd+m) + insertMatrixC;
                        Gnew(:,nV+stateInd+1:nV+stateInd+m) = Gnew(:,nV+stateInd+1:nV+stateInd+m) + insertMatrixG;
                        stateInd = stateInd + m;
                    end
                    Cnew(:,excludeStateInd+nV)=[];
                    Cnew(excludeStateInd+nV,:)=[];
                    Gnew(:,excludeStateInd+nV)=[];
                    Gnew(excludeStateInd+nV,:)=[];
                    Bnew(excludeStateInd+nV,:)=[];
                end
                
                Anew = Gnew + (Cnew - Cprev)/h;
                
                %Differential Lyapunov equation matrices, i.e., dK = E*K + K*E' + F*F'
                E = - Cnew\Anew;
                F = - Cnew\Bnew;
                Q = F*F';
                
                %Backward-Euler method applied to the Lyapunov equation
                if strcmp(thisSolver.solverType, 'TRPZ-BE')
                    
                    if strcmp(thisSolver.stateSpaceType, 'standard')
                        Pr = E - eye(nV+nStates)/(2*h);
                    elseif strcmp(thisSolver.stateSpaceType, 'reduced')
                        Pr = E - eye(nV+nStates-nTr)/(2*h);
                    end
                    Qr = Q + thisSolver.Kprev/h;
                    thisSolver.Kcur = lyap(Pr, Qr);
                    
                    Kdiag = diag(thisSolver.Kcur);
                    thisSolver.K(:,i) = Kdiag(Kind);
                    %Runge-Kutta method applied to the Lyapunov equation
                elseif strcmp(thisSolver.solverType, 'TRPZ-RK')
                    
                    if strcmp(thisSolver.stateSpaceType, 'standard')
                        
                        R1 = thisSolver.Kprev*E' + Q;
                        
                        S1 = (eye(nV+nStates) + h*E)*thisSolver.Kprev + h*R1;
                        
                        R2 = 0.5*(S1 + thisSolver.Kprev)*E' + Q;
                        
                        thisSolver.Kcur = (eye(nV+nStates) + h*E + (h^2/2)*(E^2))*thisSolver.Kprev + (h*eye(nV+nStates) + (h^2/2)*E)*R2;
                        
                        Kdiag = diag(thisSolver.Kcur);
                        thisSolver.K(:,i) = Kdiag(Kind);
                        
                    elseif strcmp(thisSolver.stateSpaceType, 'reduced')
                        
                        R1 = thisSolver.Kprev*E' + Q;
                        
                        S1 = (eye(nV+nStates-nTr) + h*E)*thisSolver.Kprev + h*R1;
                        
                        R2 = 0.5*(S1 + thisSolver.Kprev)*E' + Q;
                        
                        thisSolver.Kcur = (eye(nV+nStates-nTr) + h*E + (h^2/2)*(E^2))*thisSolver.Kprev + (h*eye(nV+nStates-nTr) + (h^2/2)*E)*R2;
                        
                        Kdiag = diag(thisSolver.Kcur);
                        thisSolver.K(:,i) = Kdiag(Kind);
                    end
                    
                end
                
                if floor(t/h) == floor(tref/h)
                    thisSolver.Krefcur = thisSolver.Kcur;
                    trefInd = i;
                    
                    Krefdiag = diag(thisSolver.Krefcur);
                    thisSolver.Kref(:,i-trefInd+1) = Krefdiag(Kind);
                elseif t > tref
                    if strcmp(thisSolver.stateSpaceType, 'standard')
                        thisSolver.Krefcur = (-(E - eye(nV+nStates)/h)\thisSolver.Krefprev'/h)';
                    elseif strcmp(thisSolver.stateSpaceType, 'reduced')
                        thisSolver.Krefcur = (-(E - eye(nV+nStates-nTr)/h)\thisSolver.Krefprev'/h)';
                    end
                    
                    Krefdiag = diag(thisSolver.Krefcur);
                    thisSolver.Kref(:,i-trefInd+1) = Krefdiag(Kind);
                end
                
                
                thisSolver.Kprev = thisSolver.Kcur;
                thisSolver.Krefprev = thisSolver.Krefcur;
                Cprev = Cnew;
                
                % save solver variables
                thisSolver.ynow = v;
                thisSolver.Y(:,i) = v;
                thisSolver.S(:,i) = s;
                thisSolver.PROPENSITYVECTOR(:,i) = propensityVector;
                thisSolver.numNewtonIters(i) = numIter;
                if calcCurrents == true
                    [SVt, SCt] = circ.getSourceVariables(v, thisSolver.tnow,...
                        QprevFull, thisSolver.T(i-1));
                    thisSolver.SC(:,i) = SCt;
                    thisSolver.SV(:,i) = SVt;
                end
                
                xprdct = [v; s];
                
                if floor(i/numT*100) > floor((i-1)/numT*100)
                    thisSolver.displayTime('override');
                end
            end
            
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %                    User Display functions                       %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function displaySolution(thisSolver)
            %plot results
            nV = thisSolver.numVars;
            nVV = thisSolver.numVoltageVars;
            figW = 562;
            figH = 504;
            scrnSize = get(0,'ScreenSize');
            scrnH = scrnSize(4);
            figPos = [150, scrnH-figH-50, figW, figH];
            
            figure('OuterPosition', figPos);
            plot(thisSolver.T,thisSolver.Y(1:nVV,:)); grid on;
            title('Noiseless membrane voltages');
            
            %switch thisSolver.currentVariables
            %    case 'on'
            %        figure('OuterPosition', figPos + [figW, 0, 0, 0]);
            %        plot(thisSolver.T,thisSolver.Y(nVV+1:nV,:)); grid on;
            %        title('Currents');
            %    case 'off'
            %        figure('OuterPosition', figPos + [figW, 0, 0, 0]);
            %        plot(thisSolver.T,thisSolver.SC); grid on;
            %        title('Source Currents');
            %end
            
            if thisSolver.numMCs > 0
                figure('OuterPosition', figPos + [figW, 0, 0, 0]);
                plot(thisSolver.T,thisSolver.S); grid on;
                title('Noiseless ion channel states');
                %if ~strncmp(thisSolver.solverType, 'IDA', 3)
                %    figure('OuterPosition', figPos + [figW, -figH, 0, 0]);
                %    plot(thisSolver.T,thisSolver.numNewtonIters); grid on;
                %    title('Newton iterations');
                %end
            end
        end
        
        function displayTime(thisSolver,varargin)
            if nargin > 1 && strcmp(varargin{1}, 'override')
                fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b\b\b%12e\n', thisSolver.tnow);
            else
                fprintf(1,'time: %12e\n', thisSolver.tnow);
            end
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %   Residual and Jacobian functions for various solvers           %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [R, flag] = IDAResFunEMS(thisSolver, t, y, yp)
            %NOTE: hard-code this scaling stuff later. You are recomputing
            %the scaling factors at every iteration because at the testing
            %stage you want to be able to change stuff without reloading
            %the classes!!!!! Do the hard-coding with sparse matrices.
            nV = thisSolver.numVars;
            nMCs = thisSolver.numMCs;
            nStates = thisSolver.numStates;
            
            Kchrg = thisSolver.chargeScaleFactor;
            
            v = y(1:nV);
            q = y(nV+1:2*nV);
            qp = yp(nV+1:2*nV);
            
            s = y(2*nV+1:2*nV+nStates);
            sp = yp(2*nV+1:2*nV+nStates);
            
            
            [Q, I, J] = thisSolver.stampFun(v,t,s);
            [rateVector, dRateVector, dRateVectorState, ~, ~, stateChangeMatrix] = thisSolver.ttcFun(v,t,s);
            
            [propensityVector] = thisSolver.propensityFun(rateVector, dRateVector, dRateVectorState, s);
            
            R1 = (qp/Kchrg + I - J);
            R2 = q - Q*Kchrg;
            R3 = sp - stateChangeMatrix*propensityVector;
            
            k=0;
            for i=1:1:nMCs
                m = thisSolver.circuit.MCs{i}.numStates;
                R3(k+m) = sum(s(k+1:k+m)) - sum(thisSolver.circuit.MCs{i}.stateVector);
                k=k+m;
            end
            
            R = [R1;
                R2;
                R3];
            flag = 0;
        end
        
        function [J, flag] = IDADJacFunEMS(thisSolver, t, y, ~, ~, cj)
            nV = thisSolver.numVars;
            nMCs = thisSolver.numMCs;
            nTransitions = thisSolver.numTransitions;
            nStates = thisSolver.numStates;
            
            Kchrg = thisSolver.chargeScaleFactor;
            
            v = y(1:nV);
            s = y(2*nV+1:2*nV+nStates);
            
            [~, ~, ~, dQ, dI] = thisSolver.stampFun(v,t,s);
            
            [rateVector, dRateVector, dRateVectorState, dQtr, dItr, stateChangeMatrix] = ...
                thisSolver.ttcFun(v,t,s);
            
            [~, dPropensityVector, dPropensityVectorStatePart1, dPropensityVectorStatePart2] = ...
                thisSolver.propensityFun(rateVector, dRateVector, dRateVectorState, s);
            
            dPropensityVectorState = dPropensityVectorStatePart1 + dPropensityVectorStatePart2;
            
            J11 = dI;
            J12 = sparse(nV,nV);
            J13 = dItr;
            
            J21 = -dQ*Kchrg;
            J22 = speye(nV);
            J23 = -dQtr*Kchrg;
            
            J31 = -stateChangeMatrix*dPropensityVector;
            J32 = sparse(nStates,nV);
            J33 = -stateChangeMatrix*dPropensityVectorState;
            
            k=0;
            H = speye(nStates,nStates);
            for i=1:1:nMCs
                m = thisSolver.circuit.MCs{i}.numStates;
                J31(k+m,:) = zeros(nV,1);
                J33(k+m,:) = zeros(nStates,1);
                J33(k+m,k+1:k+m) = ones(m,1);
                H(k+m,k+m) = 0;
                k=k+m;
            end
            
            J = full([J11, J12, J13;
                J21, J22, J23;
                J31, J32, J33] + cj*[sparse(nV,nV),  speye(nV)/Kchrg,  sparse(nV,nStates);
                sparse(nV,nV),  sparse(nV,nV),    sparse(nV,nStates);
                sparse(nStates,nV), sparse(nStates,nV),   H]);
            
            flag = 0;
            j=0;
        end
        
        function [Jv, flag] = IDAJacTimesVecFunEMS(thisSolver, t, y, ~, ~, vec, cj)
            nV = thisSolver.numVars;
            nMCs = thisSolver.numMCs;
            nTransitions = thisSolver.numTransitions;
            nStates = thisSolver.numStates;
            
            Kchrg = thisSolver.chargeScaleFactor;
            
            v = y(1:nV);
            %q = y(nV+1:2*nV);
            s = y(2*nV+1:2*nV+nStates);
            
            [~, ~, ~, dQ, dI] = thisSolver.stampFun(v,t,s);
            
            [rateVector, dRateVector, dRateVectorState, dQtr, dItr, stateChangeMatrix] = ...
                thisSolver.ttcFun(v,t,s);
            
            [~, dPropensityVector, dPropensityVectorStatePart1, dPropensityVectorStatePart2] = ...
                thisSolver.propensityFun(rateVector, dRateVector, dRateVectorState, s);
            
            dPropensityVectorState = dPropensityVectorStatePart1 + dPropensityVectorStatePart2;
            
            J11 = dI;
            J12 = sparse(nV,nV);
            J13 = dItr;
            
            J21 = -dQ*Kchrg;
            J22 = speye(nV);
            J23 = -dQtr*Kchrg;
            
            J31 = -stateChangeMatrix*dPropensityVector;
            J32 = sparse(nStates,nV);
            J33 = -stateChangeMatrix*dPropensityVectorState;
            
            k=0;
            H = speye(nStates,nStates);
            for i=1:1:nMCs
                m = thisSolver.circuit.MCs{i}.numStates;
                J31(k+m,:) = zeros(nV,1);
                J33(k+m,:) = zeros(nStates,1);
                J33(k+m,k+1:k+m) = ones(m,1);
                H(k+m,k+m) = 0;
                k=k+m;
            end
            
            J = [J11, J12, J13;
                J21, J22, J23;
                J31, J32, J33] + cj*[sparse(nV,nV),  speye(nV)/Kchrg,  sparse(nV,nStates);
                sparse(nV,nV),  sparse(nV,nV),    sparse(nV,nStates);
                sparse(nStates,nV), sparse(nStates,nV),   H];
            
            Jv = J*vec;
            flag = 0;
        end
        
        function [z, flag] = IDAPrecSolveFunEMS(thisSolver, t, y, ~, ~, r, cj)
            nV = thisSolver.numVars;
            nMCs = thisSolver.numMCs;
            nTransitions = thisSolver.numTransitions;
            nStates = thisSolver.numStates;
            
            Kchrg = thisSolver.chargeScaleFactor;
            
            v = y(1:nV);
            %q = y(nV+1:2*nV);
            s = y(2*nV+1:2*nV+nStates);
            
            [~, ~, ~, dQ, dI] = thisSolver.stampFun(v,t,s);
            
            [rateVector, dRateVector, dRateVectorState, dQtr, dItr, stateChangeMatrix] = ...
                thisSolver.ttcFun(v,t,s);
            
            [~, dPropensityVector, dPropensityVectorStatePart1, dPropensityVectorStatePart2] = ...
                thisSolver.propensityFun(rateVector, dRateVector, dRateVectorState, s);
            
            dPropensityVectorState = dPropensityVectorStatePart1 + dPropensityVectorStatePart2;
            
            
            J11 = dI;
            J12 = sparse(nV,nV);
            J13 = dItr;
            
            J21 = -dQ*Kchrg;
            J22 = speye(nV);
            J23 = -dQtr*Kchrg;
            
            J31 = -stateChangeMatrix*dPropensityVector;
            J32 = sparse(nStates,nV);
            J33 = -stateChangeMatrix*dPropensityVectorState;
            
            k=0;
            H = speye(nStates,nStates);
            for i=1:1:nMCs
                m = thisSolver.circuit.MCs{i}.numStates;
                J31(k+m,:) = zeros(nV,1);
                J33(k+m,:) = zeros(nStates,1);
                J33(k+m,k+1:k+m) = ones(m,1);
                H(k+m,k+m) = 0;
                k=k+m;
            end
            
            J = [J11, J12, J13;
                J21, J22, J23;
                J31, J32, J33] + cj*[sparse(nV,nV),  speye(nV)/Kchrg,  sparse(nV,nStates);
                sparse(nV,nV),  sparse(nV,nV),    sparse(nV,nStates);
                sparse(nStates,nV), sparse(nStates,nV),   H];
            
            z = J\r;
            flag = 0;
        end
        
        
        
        function resFun = TRPZResFun(thisSolver, t, Qprev, Iprev, Jprev, h, sprev, propensityVectorPrev)
            nV = thisSolver.numVars;
            nTr = thisSolver.numMCs;
            nStates = thisSolver.numStates;
            nTransitions = thisSolver.numTransitions;
            stateChangeMatrix = thisSolver.circuit.stateChangeMatrixi;
            stateChangeMatrixDiff = thisSolver.circuit.stateChangeMatrixiDiff;
            Kchrg = thisSolver.chargeScaleFactor;
            
            function R = res(x)
                % v... new voltages and possibly currents
                % s... new trap states
                v = x(1:nV);
                s = x(nV+1:nV+nStates);
                %seq = x(nV+nStates+1:nV+2*nStates);
                
                %sprev = sallprev(1:nStates);
                %seqprev = sallprev(nStates+1:2*nStates);
                
                [Q, I, J] = thisSolver.stampFun(v,t,s);
                %[propensityVector] = thisSolver.ttcFun(v,t,s);
                %[propensityVectorEq] = thisSolver.ttcFun(v,t,seq);
                
                [rateVector, dRateVector, dRateVectorState] = ...
                    thisSolver.ttcFun(v,t,s);
                
                [propensityVector] = ...
                    thisSolver.propensityFun(rateVector, dRateVector, dRateVectorState, s);
                %[propensityVectorEq] = ...
                %thisSolver.propensityFun(rateVector, dRateVector, dRateVectorState, seq);
                %[propensityVectorEqNow] = ...
                %thisSolver.propensityFun(rateVector, dRateVector, dRateVectorState, seqprev);
                
                %DW = stateChangeMatrixDiff*(sqrt(propensityVectorEqPrev).*dW);
                
                R1 = Q - Qprev - h/2*(J - I + Jprev - Iprev);
                
                R2 = s - sprev - ...
                    (h/2)*(stateChangeMatrix*(propensityVector + propensityVectorPrev)); % + DW;
                
                R = [Kchrg*R1; R2];
            end
            
            resFun = @res;
        end
        
        function jacFun = TRPZJacFun(thisSolver, t, h)
            nV = thisSolver.numVars;
            nStates = thisSolver.numStates;
            nTransitions = thisSolver.numTransitions;
            stateChangeMatrix = thisSolver.circuit.stateChangeMatrixi;
            Kchrg = thisSolver.chargeScaleFactor;
            
            function J = jac(x)
                % v... new voltages and possibly currents
                % s... new trap states
                v = x(1:nV);
                s = x(nV+1:nV+nStates);
                %seq = x(nV+nStates+1:nV+2*nStates);
                
                [~, ~, ~, dQ, dI] = thisSolver.stampFun(v,t,s);
                %[~, dPropensityVector, dPropensityVectorState, dQst, dIst] = ...
                %    thisSolver.ttcFun(v,t,s);
                
                %[~, dPropensityVectorEq, dPropensityVectorStateEq] = ...
                %    thisSolver.ttcFun(v,t,seq);
                
                [rateVector, dRateVector, dRateVectorState, dQst, dIst] = ...
                    thisSolver.ttcFun(v,t,s);
                
                [~, dPropensityVector, dPropensityVectorStatePart1, dPropensityVectorStatePart2] = ...
                    thisSolver.propensityFun(rateVector, dRateVector, dRateVectorState, s);
                dPropensityVectorState =  dPropensityVectorStatePart1 + dPropensityVectorStatePart2;
                
                %[~, dPropensityVectorEq, dPropensityVectorStateEqPart1, dPropensityVectorStateEqPart2] = ...
                %thisSolver.propensityFun(rateVector, dRateVector, dRateVectorState, seq);
                
                J11 = dQ + h/2*dI;
                J12 = dQst + h/2*dIst;
                
                J21 = -(h/2)*(stateChangeMatrix*dPropensityVector);
                J22 = eye(nStates) - (h/2)*(stateChangeMatrix*dPropensityVectorState);
                
                J = [Kchrg*J11, Kchrg*J12;
                    J21,       J22];
            end
            jacFun = @jac;
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %                    Non-linear solver                            %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [x, numIter] = newtSol(thisSolver, inFun, dInFun, x0, varargin)
            % NEWTSOL Newton's method for the solution of non-linear equations
            % inFun... function handle to the residual function
            % dInFun... Jacobian of the residual function
            % x0 initial guess
            % varargin can be 'Damped' for damped Newton's method
            atol = thisSolver.tolVec;
            rtol = thisSolver.reltol;
            restol = thisSolver.restol;
            damping = false;
            % is damping on?
            if (nargin > 4) && strcmp(varargin{1}, 'Damped')
                damping = true;
                imax = 4;
            end
            numIter = 0;
            convCheck = false;
            x = x0;
            inFunVal = inFun(x);
            dampInd = 1;
            numFunEval = 1;
            while convCheck == false
                dInFunVal = dInFun(x);
                try
                    lastwarn('');
                    deltax = - dInFunVal\inFunVal;
                    error(lastwarn);
                catch exception
                    rethrow(exception);
                end
                x_new = x + deltax;
                lsVecPrev = inFun(x_new);
                lsNormPrev = norm(lsVecPrev./restol);
                % if damping is on, do a line search
                if damping == true
                    for dampInd = 1:imax
                        lsVecCur = inFun(x + deltax/(2^dampInd));
                        lsNormCur = norm(lsVecCur./restol);
                        if lsNormCur > lsNormPrev
                            break;
                        end
                        lsVecPrev = lsVecCur;
                        lsNormPrev = lsNormCur;
                    end
                    alpha = 1/(2^(dampInd-1));
                    x_new = x + alpha*deltax;
                end
                numFunEval = numFunEval + dampInd;
                inFunVal = lsVecPrev;
                convCheck = (norm((x_new - x)./(atol + rtol*abs(x_new))) <= 1) && ...
                    (lsNormPrev <= 1);
                x = x_new;
                numIter = numIter + 1;
            end
        end
        
        
        function [propensityVector, dPropensityVector, dPropensityVectorStatePart1, dPropensityVectorStatePart2] ...
                = propensityFun(thisSolver, rateVector, dRateVector, dRateVectorState, s)
            
            m=0;
            n=0;
            k=m;
            l=n;
            thisCircuitMC = thisSolver.circuit;
            numCircuitStates = thisCircuitMC.numStates;
            nMCs = thisCircuitMC.numMCs;
            parentComponentCurrent = thisCircuitMC.MCs{1}.parentComponent;
            statesCircComp = parentComponentCurrent.numStates;
            transitionsCircComp = parentComponentCurrent.numTransitions;
            for j=1:nMCs
                
                MC = thisCircuitMC.MCs{j};
                
                if isequal(parentComponentCurrent, MC.parentComponent) == false
                    parentComponentPrev = parentComponentCurrent;
                    parentComponentCurrent = MC.parentComponent;
                    m = m + parentComponentPrev.numTransitions;
                    n = n + parentComponentPrev.numStates;
                    k = m;
                    l = n;
                end
                
                statesMC = MC.numStates;
                transitionsMC = MC.numTransitions;
                
                [sPropensityVector, sdPropensityVector, sdPropensityVectorStatePart1, sdPropensityVectorStatePart2] ...
                    = MC.propensityEval(rateVector(k+1:k+transitionsMC,1), dRateVector(k+1:k+transitionsMC,:), dRateVectorState(k+1:k+transitionsMC,:), s(l+1:l+statesMC,1), l, numCircuitStates);
                
                propensityVector(k+1:k+transitionsMC,1) = sPropensityVector;
                dPropensityVector(k+1:k+transitionsMC,:) = sdPropensityVector;
                dPropensityVectorStatePart1(k+1:k+transitionsMC,:) = sdPropensityVectorStatePart1;
                dPropensityVectorStatePart2(k+1:k+transitionsMC,:) = sdPropensityVectorStatePart2;
                k=k+transitionsMC;
                l=l+statesMC;
            end
        end
        
    end
    
    methods (Static)
        function [BigS, dLc, dLe] = expandTtc(s, Le, Lc, dTc, dTe, N)
            BigS = s(:,ones(1,N));
            BigLe = Le(:,ones(1,N));
            BigLc = Lc(:,ones(1,N));
            dLc = -((BigLc).^2).*dTc;
            dLe = -((BigLe).^2).*dTe;
        end
    end
    
    % end methods
end