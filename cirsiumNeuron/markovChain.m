classdef markovChain < handle
    % MARKOVCHAIN Abstract class for Markov Chain models.

    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------
    
    properties (SetAccess = protected)
        stateVector;        % state vector
        numStates;          % number of states
        numTransitions;     % number of transitions
        stateChangeMatrix;  % state change vectors in matrix form
        rateVector;         % state change rate vector
        dRateVector;
        dRateVectorState;
        dQst;
        dIst;
        lambda;             % total rate for the next event
        dLambda;            % derivative of lambda wrt voltages
        number;             % number of the Markov Chain in the parent circuit
    end
    
    properties
        updateNeeded = true;
    end
    
    
    methods (Abstract)
        setState(thisMarkovChain, newStateVector);
        flipState(thisMarkovChain, indEvent);
        updateRates(thisMarkovChain);
        rateVector = getRates(thisMarkovChain);
    end
    
 
    methods
        function thisMarkovChain = markovChain(stateVector,stateChangeMatrix)
            % constructor method for a Markov Chain.
            thisMarkovChain.stateVector = stateVector;
            thisMarkovChain.stateChangeMatrix = stateChangeMatrix;  
        end
        
        function setNumber(thisMarkovChain, n)
            thisMarkovChain.number = n;
        end
    end
    
end

