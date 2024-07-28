classdef NaChannelMC < markovChain
    % NACHANNELMC creates Na channels to be embedded in neuronMC.

    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------

    properties (SetAccess = private)
        
        g_Na = 20*10^-12; %single Na channel conductance
        
        parentComponent = resistor.empty;       
        
    end
    
    properties
        eqNums; % nums of top level variables connected to the trap
        neighborNums; % top level equation nums of neighbor traps
    end
    
    methods
        function thisNaChannel = NaChannelMC(stateVector, stateChangeMatrix, parentNeuron)
            
            thisNaChannel = thisNaChannel@markovChain(stateVector, stateChangeMatrix);
            
            thisNaChannel.parentComponent = parentNeuron;
            
            thisNaChannel.numStates = 8; 
            thisNaChannel.numTransitions = 20;
            
        end
        
        

        function setState(thisNaChannel, newStateVector)
            
            oldStateVector = thisNaChannel.stateVector; 

            if any(newStateVector ~= oldStateVector)
                thisNaChannel.stateVector = newStateVector;
                thisNaChannel.parentComponent.updateMCEffects(thisNaChannel, oldStateVector, newStateVector);
                thisNaChannel.parentComponent.forceEquationGeneration();
            end  
        end
        
        
        
        function flipState(thisNaChannel, indEvent)
            
            oldStateVector = thisNaChannel.stateVector;

            newStateVector = thisNaChannel.stateVector + thisNaChannel.stateChangeMatrix(:,indEvent);
            
            if sum(oldStateVector)~=sum(newStateVector)
               
                error('The number of traps are not the same in the old and new state vectors!'); 
            end 
            
            setState(thisNaChannel, newStateVector);
              
        end
        
        
        function updateRates(thisNaChannel)
            
            [rateVector, dRateVector, dRateVectorState, dQst, dIst] = thisNaChannel.parentComponent.MCRates(thisNaChannel);

            thisNaChannel.rateVector = rateVector;
            thisNaChannel.dRateVector = dRateVector;
            thisNaChannel.dRateVectorState = dRateVectorState;
            thisNaChannel.dQst = dQst;
            thisNaChannel.dIst = dIst;
            
            thisNaChannel.updateNeeded = false;
        end

        
        function [rateVector, dRateVector, dRateVectorState, dQst, dIst] = getRates(thisNaChannel)
            
            if thisNaChannel.updateNeeded == true
                thisNaChannel.updateRates();
            end
            
            rateVector = thisNaChannel.rateVector;
            dRateVector = thisNaChannel.dRateVector;
            dRateVectorState = thisNaChannel.dRateVectorState;
            dQst = thisNaChannel.dQst;
            dIst = thisNaChannel.dIst;
        end
        
        
        function [propensityVector, dPropensityVector, dPropensityVectorStatePart1, dPropensityVectorStatePart2] ...
                = propensityEval(thisNaChannel, rateVector, dRateVector, dRateVectorState, stateVector, N, numCircuitStates)
            
            nStates = numCircuitStates;
            stateChangeMatrix = thisNaChannel.stateChangeMatrix;
            
            m = thisNaChannel.numStates;
            n = thisNaChannel.numTransitions;
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