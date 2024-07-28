classdef KChannelMC < markovChain
    % KCHANNELMC creates K channels to be embedded in neuronMC.

    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------
    
    properties (SetAccess = private)
        
        g_K  = 20*10^-12; %single K channel conductance
        parentComponent = resistor.empty;       

    end
    
    properties
        eqNums;         % nums of top level variables connected to the trap
        neighborNums;   % top level equation nums of neighbor traps
    end
    
    methods
        function thisKChannel = KChannelMC(stateVector, stateChangeMatrix, parentNeuron)
            
            thisKChannel = thisKChannel@markovChain(stateVector, stateChangeMatrix);
            
            thisKChannel.parentComponent = parentNeuron;
            
            thisKChannel.numStates = 5; 
            thisKChannel.numTransitions = 8;
            
        end
        
        

        function setState(thisKChannel, newStateVector)
            
            oldStateVector = thisKChannel.stateVector;  

            if any(newStateVector ~= oldStateVector)
                thisKChannel.stateVector = newStateVector;
                thisKChannel.parentComponent.updateMCEffects(thisKChannel, oldStateVector, newStateVector);
                thisKChannel.parentComponent.forceEquationGeneration();
            end 
        end
        
        
        
        function flipState(thisKChannel, indEvent)
            oldStateVector = thisKChannel.stateVector;

            newStateVector = thisKChannel.stateVector + thisKChannel.stateChangeMatrix(:,indEvent);
            
            if sum(oldStateVector)~=sum(newStateVector)
               
                error('The number of traps are not the same in the old and new state vectors!'); 
            end
            
            setState(thisKChannel, newStateVector);
             
        end
        
        
        function updateRates(thisKChannel)
            
            [rateVector, dRateVector, dRateVectorState, dQst, dIst] = thisKChannel.parentComponent.MCRates(thisKChannel);
            
            thisKChannel.rateVector = rateVector;
            thisKChannel.dRateVector = dRateVector;
            thisKChannel.dRateVectorState = dRateVectorState;
            thisKChannel.dQst = dQst;
            thisKChannel.dIst = dIst;
            
            thisKChannel.updateNeeded = false;
        end

        
        function [rateVector, dRateVector, dRateVectorState, dQst, dIst] = getRates(thisKChannel)
            
            if thisKChannel.updateNeeded == true
                thisKChannel.updateRates();
            end
            
            rateVector = thisKChannel.rateVector;
            dRateVector = thisKChannel.dRateVector;
            dRateVectorState = thisKChannel.dRateVectorState;
            dQst = thisKChannel.dQst;
            dIst = thisKChannel.dIst;
        end
        
        function [propensityVector, dPropensityVector, dPropensityVectorStatePart1, dPropensityVectorStatePart2] ...
                = propensityEval(thisKChannel, rateVector, dRateVector, dRateVectorState, stateVector, N, numCircuitStates)
            
            nStates = numCircuitStates;
            stateChangeMatrix = thisKChannel.stateChangeMatrix;
            
            m = thisKChannel.numStates;
            n = thisKChannel.numTransitions;
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