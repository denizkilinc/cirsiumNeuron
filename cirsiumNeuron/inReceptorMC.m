classdef inReceptorMC < markovChain
    % INRECEPTORMC creates inhibitory receptors to be embedded in neuronMC.

    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------

    properties (SetAccess = private)
        
        dRateVectorCross;
        g_InReceptor = 20*10^-12; %single receptor channel conductance
        
        parentComponent = resistor.empty;
        
    end
    
    properties
        sourceNeuron;
        eqNums;         % nums of top level variables connected to the trap
        neighborNums;   % top level equation nums of neighbor traps
    end
    
    methods
        function thisReceptor = inReceptorMC(stateVector, stateChangeMatrix, parentNeuron)
            
            thisReceptor = thisReceptor@markovChain(stateVector, stateChangeMatrix);
            
            thisReceptor.parentComponent = parentNeuron;
            
            thisReceptor.numStates = 5;
            thisReceptor.numTransitions = 8;
            
        end
        
        
        function setState(thisReceptor, newStateVector)
            
            oldStateVector = thisReceptor.stateVector;
            
            if any(newStateVector ~= oldStateVector)
                thisReceptor.stateVector = newStateVector;
                thisReceptor.parentComponent.updateMCEffects(thisReceptor, oldStateVector, newStateVector);
                thisReceptor.parentComponent.forceEquationGeneration();
            end
        end
        
        
        
        function flipState(thisReceptor, indEvent)
            oldStateVector = thisReceptor.stateVector;
            
            newStateVector = thisReceptor.stateVector + thisReceptor.stateChangeMatrix(:,indEvent);
            
            if sum(oldStateVector)~=sum(newStateVector)
                
                error('The number of traps are not the same in the old and new state vectors!');
            end
            
            setState(thisReceptor, newStateVector);
        end
        
        
        function updateRates(thisReceptor)
            
            [rateVector, dRateVector, dRateVectorState, dQst, dIst] = thisReceptor.parentComponent.MCRates(thisReceptor);
            
            thisReceptor.rateVector = rateVector;
            thisReceptor.dRateVector = dRateVector;
            thisReceptor.dRateVectorState = dRateVectorState;
            thisReceptor.dQst = dQst;
            thisReceptor.dIst = dIst;
            
            sourceNeuron = thisReceptor.sourceNeuron;
            sinkNeuron = thisReceptor.parentComponent;
            if isa(sourceNeuron, 'neuronMC') == true
                [dRateVectorCross] = synapseRates(sourceNeuron, sinkNeuron, thisReceptor);
            elseif isa(sourceNeuron, 'neuronHH') == true
                [dRateVectorCross] = synapseRates(sinkNeuron, sourceNeuron, thisReceptor);
            end
            thisReceptor.dRateVectorCross = dRateVectorCross;
            
            thisReceptor.updateNeeded = false;
        end
        
        
        function [rateVector, dRateVector, dRateVectorState, dQst, dIst] = getRates(thisReceptor)
            
            if thisReceptor.updateNeeded == true
                thisReceptor.updateRates();
            end
            
            rateVector = thisReceptor.rateVector;
            dRateVector = thisReceptor.dRateVector;
            dRateVectorState = thisReceptor.dRateVectorState;
            dQst = thisReceptor.dQst;
            dIst = thisReceptor.dIst;
        end
        
        
        
        function [dRateVectorCross] = getCrossRates(thisReceptor)
            
            sourceNeuron = thisReceptor.sourceNeuron;
            
            if thisReceptor.updateNeeded == true
                thisReceptor.updateRates();
            end
            
            dRateVectorCross = thisReceptor.dRateVectorCross;
        end
        
        
        function [propensityVector, dPropensityVector, dPropensityVectorStatePart1, dPropensityVectorStatePart2] ...
                = propensityEval(thisReceptor, rateVector, dRateVector, dRateVectorState, stateVector, N, numCircuitStates)
            
            nStates = numCircuitStates;
            stateChangeMatrix = thisReceptor.stateChangeMatrix;
            
            
            m = thisReceptor.numStates;
            n = thisReceptor.numTransitions;
            
            [a, b]=size(dRateVector);
            
            propensityVector = zeros(n,1);
            dPropensityVector = zeros(a,b);
            dPropensityVectorStatePart1 = zeros(n, nStates);
            dPropensityVectorStatePart2 = zeros(n, nStates);
            
            for i=1:1:m
                
                for j=1:1:n
                    
                    if stateChangeMatrix(i,j) == -1
                        
                        propensityVector(j,1) = rateVector(j,1)*stateVector(i,1);
                        
                        dPropensityVector(j,:) = dRateVector(j,:)*stateVector(i,1);
                        
                        dPropensityVectorStatePart1(j,:) = dRateVectorState(j,:)*stateVector(i,1);
                        
                        dPropensityVectorStatePart2(j,i+N) = rateVector(j,1);
                        
                    end
                end
            end
        end
        
        
    end
end

