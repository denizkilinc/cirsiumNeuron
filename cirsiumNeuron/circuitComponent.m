classdef circuitComponent < electrical & matlab.mixin.Copyable & matlab.mixin.Heterogeneous
    % CIRCUITCOMPONENT Container and controller class for circuit
    % components.
    
    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------
    
    properties (SetAccess = protected)
        terminals = node.empty; % Terminals in form of node references.
        parentCircuit = circuit.empty; % The circuit this device belongs to.
        currentVarInd = 0;                 % The index of the current variable in the whole circuit.
        voltEqNums;
        currEqNums;
        termNodeNumbers;    % Node numbers of the terminals in the circuit.
    end
    
    methods (Abstract)
        sQ = Q(thisCircuitComponent, x, t); % Q vector of the device.
        sI = I(thisCircuitComponent, x, t); % I vector of the device.
        sJ = J(thisCircuitComponent, x, t); % J vector of the device.
        sdQ = dQ(thisCircuitComponent, x, t); % Jacobian of the Q vector.
        sdI = dI(thisCircuitComponent, x, t); % Jacobian of the I vector.
        setEquationNumbers(thisCircuitComponent);
        tlc = getTopLevelCircuit(thisCircuitComponent);
    end
    
    methods 
        function obj = circuitComponent(name)
            % CIRCUITELEMENT constructor for a circuit component.
            if nargin > 0
                obj.name = name;
            end
        end
        
        function sealComponent(thisCircuitComponent)
            if isempty(thisCircuitComponent.parentCircuit)
                error('This circuit component has no parent circuit!');
            end
            
            for i = 1:thisCircuitComponent.numComponents
                cc = thisCircuitComponent.components(i);
                cc.sealComponent;
                if cc.numVoltageSources > 0
                    thisCircuitComponent.appendVoltageSource(cc);
                end
            end
            thisCircuitComponent.termNodeNumbers = ...
                                    [thisCircuitComponent.terminals.number];
            thisCircuitComponent.setEquationNumbers();
            thisCircuitComponent.isSealed = true;
        end
        
        
        function setParent(thisCircuitComponent, circ)
            if ~isempty(thisCircuitComponent.parentCircuit)
                error('This device already has a parent!');
            end
            
            thisCircuitComponent.parentCircuit = circ;            
            % If the device needs additional current variables such as a
            % voltage source or an inductor
            if thisCircuitComponent.numCurrentVars ~= 0
                thisCircuitComponent.currentVarInd = circ.numCurrentVars + 1;
            end
        end

        function setTerminal(thisCircuitComponent, termNum, aNode)
            N = thisCircuitComponent.numTerms;
            if termNum < 0 || termNum > N
                error('No such terminal number');
            elseif ~isa(aNode, 'node')
                error('Terminals must be of type "node"!');
            end                         
            
            thisCircuitComponent.terminals(termNum) = aNode;
            aNode.addComponent(thisCircuitComponent);
        end
    end
    
    methods (Access = protected)
        function newCC = copyElement(thisCircuitComponent)
            newCC = copyElement@matlab.mixin.Copyable(thisCircuitComponent);
            newCC.parentCircuit = circuit.empty;
            newCC.terminals = node.empty;
            newCC.groundNode = node.empty;
            newCC.numComponents = 0;
            newCC.numVoltageVars = 0;
            newCC.voltEqNums = [];
            newCC.currEqNums = [];
            newCC.terminals = node.empty;
            newCC.parentCircuit = circuit.empty;
            newCC.currentVarInd = 0;
            newCC.voltageSources = voltageSource.empty;
            newCC.numVoltageSources = 0;
            newCC.voltageSourceNodeNums = [];
        end
    end
    
    methods (Sealed)
        function isEq = eq(comp1, comp2)
            isEq = eq@handle(comp1,comp2);
        end
        
        function obj = findobj(comps, propName, propVal)
            obj = findobj@handle(comps, propName, propVal);
        end
    end
    
    methods (Static, Sealed, Access = protected) 
        function defObj = getDefaultScalarElement()
            defObj = resistor();
        end
    end       

end