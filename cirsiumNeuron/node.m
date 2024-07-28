classdef node < handle
    % NODE A simple node class.
    % This class is a simple implementation of circuit nodes without much
    % functionality. It only keeps a list of circuit components and
    % connects them to its parent circuit.
    
    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------
        
    properties (SetAccess = private)
        name;           % Name of the node.
        numComponents = 0; % Number of devices connected to this node.
    end
    
    properties
        number = -1;    % Number of the node in its parentcircuit.
        components = resistor.empty;   % Cell array of connected devices to this node.
        parentCircuit;  % Parent circuit of the node.
    end

    methods
        function obj = node(name, number)
            %NODE constructor for node class
            if nargin > 0
                obj.name = name;
                if nargin > 1
                    obj.number = number;
                end
            end
        end
        
        function addComponent(thisNode, circComp)
            % ADDELEMENT Add a circuit component to this node.
            if ~any(thisNode.components == circComp)
                nComp = thisNode.numComponents + 1;
                thisNode.numComponents = nComp;
                thisNode.components(nComp) = circComp;
            end
        end
        
        function tln = getTopLevelNumber(thisNode)
            if ~isempty(thisNode.parentCircuit)
                tln = thisNode.parentCircuit.getTopLevelNodeNumber(thisNode.number);
            else
                error('This node does not have a parent circuit!');
            end
        end
        
        function tlc = getTopLevelCircuit(thisNode)
            if ~isempty(thisNode.parentCircuit)
                tlc = thisNode.parentCircuit.getTopLevelCircuit();
            else
                error('This node does not have a parent circuit!');
            end
        end
    end
    
    
end

