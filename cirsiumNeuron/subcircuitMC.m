classdef subcircuitMC < circuitMC & subcircuit
    % SUBCIRCUITMC Class for a circuitMC that can be added to another circuitMC.
    % This class is derived from both CIRCUTMC and SUBCIRCUIT classes.
    
    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------
    
    properties (SetAccess = protected)
        MCInd = 0;
        MCEqNums;
    end
    
    methods
        function newSubcircuitMC = subcircuitMC(varargin)
            newSubcircuitMC = newSubcircuitMC@circuitMC(varargin{:});
            newSubcircuitMC = newSubcircuitMC@subcircuit();
        end
        
        function setParent(thisSubcircuitMC, circ)
            setParent@subcircuit(thisSubcircuitMC, circ);
            % set the MC index of the subcircuit
            thisSubcircuitMC.MCInd = circ.numMCs + 1;
        end
        
        function setEquationNumbers(thisSubcircuitMC)
            setEquationNumbers@subcircuit(thisSubcircuitMC);
            MC_Ind = thisSubcircuitMC.MCInd;
            nMC = thisSubcircuitMC.numMCs;
            MCNums = MC_Ind:MC_Ind+nMC-1;
            thisSubcircuitMC.MCEqNums = MCNums;
        end
        
    end
    
    methods (Access = protected)
        function newScMC = copyElement(thisSubcircuitMC)
            nMCTemp = thisSubcircuitMC.numMCs;
            MCTemp = thisSubcircuitMC.MCs;
            thisSubcircuitMC.numMCs = 0;
            thisSubcircuitMC.MCs = markovChain.empty;
            
            newScMC = copyElement@subcircuit(thisSubcircuitMC);
            thisSubcircuitMC.numMCs = nMCTemp;
            thisSubcircuitMC.MCs = MCTemp;
            newScMC.XMCi = [];
            newScMC.VMCi = [];
            newScMC.tMCi = -Inf;
        end
    end
end

