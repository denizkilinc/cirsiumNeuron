classdef circuitComponentMC < subcircuitMC
    % CIRCUITCOMPONENTMC Container and controller class for circuit 
    % components which include subprocesses modeled by Markov Chains.

    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------

    
    properties (SetAccess = public)
        number;
    end
    
    methods (Abstract)
        [rateVector, dRateVector, dRateVectorState, dQst, dIst] = MCRates(thisComponentMC, aMC);
        registerMCEffects(thisComponentMC, aMC, oldStateVector, newStateVector);
        numMCs = determineNumberOfMCs(thisComponentMC);
        
    end
    
    methods (Abstract, Access = protected)
        initMCs(thisComponentMC, numMCs);
    end 
    
    methods
        function thisComponentMC = circuitComponentMC(varargin)
            
            thisComponentMC = thisComponentMC@subcircuitMC(varargin{:});
        end
        
        
        function addComponent(thisComponentMC, circComp, varargin)
           if thisComponentMC.numComponents == 0
               addComponent@subcircuitMC(thisComponentMC, circComp, varargin{:});
               thisComponentMC.setExternalNodes(varargin{:});
           else
               error(['This circuit component with MCs already has '...
                   'a nucleus component!']);
           end
        end
        
        function updateMCEffects(thisComponentMC, aMC, oldState, newState)
            thisComponentMC.registerMCEffects(aMC, oldState, newState);
            thisComponentMC.notifyMCs();
        end
        
    end
    
    methods (Access = protected)        
        function [Q I J dQ dI] = eqEval(thisComponentMC, x, t)
            [Q I J dQ dI] = eqEval@subcircuitMC(thisComponentMC, x, t);
            thisComponentMC.notifyMCs();
        end
        
        function notifyMCs(thisComponentMC)
            for i=1:thisComponentMC.numMCs
                thisComponentMC.MCs{i}.updateNeeded = true;
            end
        end
    end
end

