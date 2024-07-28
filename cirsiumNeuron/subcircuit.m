classdef subcircuit < circuit & circuitComponent
    % SUBCIRCUIT Class for a circuit that can be added to another circuit.
    % This class is derived from both CIRCUT and CIRCUITCOMPONENT classes.
    
    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------

    properties (SetAccess = protected)
        numTerms = 0;
        externalNodes = node.empty;     % Array of external nodes
        exNodeNumbers = 0;
        voltageVarInd = 0;
    end      
    
    methods
        function obj = subcircuit(varargin)
            obj = obj@circuit(varargin{:});
        end
        
        function sealComponent(thisSubcircuit)
            % we want to generate full equations - no ground removal!
            thisSubcircuit.groundNodeNumber = 0;
            
            % set number of independent voltage variables
            thisSubcircuit.numIndepVoltVars = thisSubcircuit.numVoltageVars -...
                                                ~isempty(thisSubcircuit.groundNode);
            
            % set total number of variables
            thisSubcircuit.numVars = thisSubcircuit.numIndepVoltVars + ...
                                                thisSubcircuit.numCurrentVars;
                                            
            % set external node numbers
            thisSubcircuit.exNodeNumbers = [thisSubcircuit.externalNodes.number];
            if isempty(thisSubcircuit.exNodeNumbers)
                error('This subcircuit does not have any external nodes!');
            end
            
            thisSubcircuit.Xi = NaN(thisSubcircuit.numVars,1);
            thisSubcircuit.Vi = NaN(thisSubcircuit.numVarsnc,1);
            
            % this has to come last
            sealComponent@circuitComponent(thisSubcircuit);
        end
    
        function setEquationNumbers(thisSubcircuit)             
                termNums = thisSubcircuit.termNodeNumbers;
                currInd = thisSubcircuit.currentVarInd;
                voltInd = thisSubcircuit.parentCircuit.numNodes + ...
                                            thisSubcircuit.voltageVarInd;
                
                
                nVV = thisSubcircuit.numVoltageVars;
                nCV = thisSubcircuit.numCurrentVars;
                nExN = thisSubcircuit.numTerms;
                
                exNodeNums = thisSubcircuit.exNodeNumbers;

                indMap1 = [exNodeNums setdiff(1:nVV, exNodeNums)];
                indMap2 = [termNums voltInd:(voltInd+nVV-1-nExN)];

                vNums = zeros(1,nVV);
                vNums(indMap1) = indMap2;

                cNums = currInd:(currInd+nCV-1);
                
                

                if ~isempty(thisSubcircuit.groundNode)
                    gGrNum = thisSubcircuit.terminals(...
                             thisSubcircuit.externalNodes ==...
                             thisSubcircuit.groundNode).number;
                    vNums(vNums == gGrNum) = [];
                end
                thisSubcircuit.voltEqNums = vNums;
                thisSubcircuit.currEqNums = cNums;       
                
        end
        
        function setExternalNodes(thisSubcircuit, varargin)
            % SETEXTERNALNODES Assign external nodes to this circuit.
            %   Subcircuits are connected to other circuits trough their
            %   external nodes. This method sets the nodes in varargin as
            %   external nodes. The nodes in varargin can be given either
            %   by their reference or name.
            N = nargin-1;
            if(N > 0)
                thisSubcircuit.externalNodes(N) = node;
                for i = 1:N
                    exNode = thisSubcircuit.findNode(varargin{i});
                    if ~isempty(exNode)
                        thisSubcircuit.externalNodes(i) = exNode;
                    else
                        error(['This node cannot be set as an external '...
                               'node because it either does not exist '...
                               'or is not part of this subcircuit!']);
                    end
                end
                thisSubcircuit.numTerms = N;
            end
        end
        
        function tlnn = getTopLevelNodeNumber(thisSubcircuit, aNodeNumber)
            % call the method recursively with the equation number this
            % subcircuit is assigned by its parent circuit.
            tlnn = thisSubcircuit.parentCircuit.getTopLevelNodeNumber(...
                    thisSubcircuit.voltEqNums(aNodeNumber));
        end
        
        function setParent(thisSubcircuit, circ)
            setParent@circuitComponent(thisSubcircuit, circ);
            
            % set the voltage index of the subcircuit
            thisSubcircuit.voltageVarInd = circ.numVoltageVars - ...
                                                        circ.numNodes + 1;
            
            % Check if the ground node of the circuit is connected to the
            % global ground
            subGround = thisSubcircuit.groundNode;
            if ~isempty(subGround)
                 exGrInd = find(thisSubcircuit.externalNodes == subGround);
                 if isempty(circ.groundNode)
                     error('The parent circuit does not have a ground node!');
                 elseif isempty(exGrInd)
                     error('Ground node must be an external node!');
                 elseif thisSubcircuit.terminals(exGrInd) ~= circ.groundNode
                     error(['Ground node must be connected to the'...
                            ' ground of the parent circuit!']);
                 end
            end
        end
        
        function tlc = getTopLevelCircuit(thisSubcircuit)
            if isempty(thisSubcircuit.parentCircuit)
                tlc = thisSubcircuit;
            else
                tlc = thisSubcircuit.parentCircuit.getTopLevelCircuit();
            end
        end     
        
        function forceEquationGeneration(thisSubcircuit)
            thisSubcircuit.forceEqGen = true;
            if isempty(thisSubcircuit.parentCircuit) == false
                thisSubcircuit.parentCircuit.forceEquationGeneration();
            end
        end
    end
    
    methods (Access = protected)
        function newSc = copyElement(thisSubcircuit)
            newSc = copyElement@circuitComponent(thisSubcircuit);
            
            newSc.nodes = node.empty;
            newSc.numNodes = 0;
            newSc.components = resistor.empty;
            newSc.numCurrentVars = 0;
            newSc.Xi = [];
            newSc.Vi = [];
            newSc.ti = -Inf;
            newSc.forceEqGen = false;
            newSc.voltageVarInd = 0;
            newSc.externalNodes = node.empty;
            newSc.voltageSourceSigns = [];
            
            for i = 1:thisSubcircuit.numNodes
                newSc.addNode(thisSubcircuit.nodes(i).name);
            end

            for i = 1:thisSubcircuit.numComponents
                cc = thisSubcircuit.components(i);
                newSc.addComponent(cc.copy, cc.terminals.name);
            end
            
            grnd = thisSubcircuit.groundNode;
            if ~isempty(grnd)
                newSc.setGroundNode(grnd.name);
            end
            newSc.setExternalNodes(thisSubcircuit.externalNodes.name);
        end
       
    end
    
end

