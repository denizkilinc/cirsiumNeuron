classdef circuit < electrical
    % CIRCUIT Container and controller class for neuronal circuits and components.
    %
    % Circuit acts as a container and controller class for neuronal circuits,
    % ion channels and synapses.
    % Circuit sets up the necessary vectors and matrices for the numerical
    % solution of the KCL and KVL equations in the following form.
    %
    %               dQ(t)/dt + I(t) = J(t)                  (1)
    %
    % Here, Q, I, J are Nx1 arrays of node charges, node currents and
    % branch source voltages respectively.
    % In addition, the jacobians of Q and I are also provided.
    % 
    % For the following discription we define
    %
    % nn ... number of nodes in the circuit
    % nc ... number of current variables in the circuit
    % ns ... number of independent voltage sources in the circuit
    %
    % Note that nc = ns + (number of inductive branches).
    % An inductive branch is usually two nodes connected by an inductor.
    %
    % For the solution of the circuit equations, the differential equation
    % (1) can be set up in two ways.
    % 1. Currents of voltage sources are treated as variables.
    % 2. Currents and the second nodes of voltage sources are eliminated.
    % 
    % In the first case, total number of variables, N, in the equation
    % system is given as follows
    %               N_1 = nn + nc
    % In the second case 
    %               N_2 = nn - ns
    %
    % To populate the circuit, use the addComponent method.
    % Sample code for an RLC circuit is given below:
    %
    %
    %       (n1)      R   (n2)   L    (n3)  C
    %         o----\/\/\/--o--(((((((--o----|(---+
    %        +|                                  |
    %        (V)                                 |
    %        -|                         <-i      |
    %         +----------------------------------o
    %                                         (gnd)
    % circ = circuit('Series RLC Circuit');
    % r = circ.addComponent(resistor(1), 'n1', 'n2');
    % l = circ.addComponent(inductor(1), 'n2', 'n3');
    % c = circ.addComponent(capacitor(1), 'n3', 'gnd');
    % v = circ.addComponent(voltageSourceDC(1), 'n1', 'gnd');
    % circuit.setGroundNode('gnd');
    % circuit.seal;
    %
    % Please see the descriptions of the individual methods for details.
    %
    % See also TRANSIENTSOLVERMC, SUBCIRCUIT, CIRCUITCOMPONENT, DEVICE, NODE.
    
    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------
        
    properties (SetAccess = protected)
        nodes = node.empty;         % An array that stores the nodes of the circuit.
        numNodes = 0;               % Number of nodes.
        numCurrentVars = 0;         % Number of independent current variables.
        numVarsnc;                  % Number of total variables - no currents case
        numIndepVoltVars = 0;       % Number of independent variables
        voltageSourceSigns;         % direction correction for voltageSourceNodeNums
        eqNumsnc;                   % equation number array for no current variable case
        Tempi = 273.15 + 27;        % Temperature in Kelvin.
    end
    
    properties (SetAccess = protected)
        % internal electrical variables
        Qi;                  % Charge array.
        Ii;                  % Current array.
        Ji;                  % Source array.
        dQi;                 % Jacobian of the charge array.
        dIi;                 % Jacobian of the current array.
        Xi;                  % Array of voltage and current variables in the circuit.
        Qinc;                % Charge array for no-current calculations.
        Iinc;                % Current array for no-current calculations.
        Jinc;                % Source array for no-current calculations (equals ZERO for now).
        dQinc;               % Jacobian of the charge array for no-current calculations.
        dIinc;               % Jacobian of the current array for no-current calculations.
        Vi;                  % Array of voltage variables to be used in current-free equations
        ti = -Inf;             % Time point in the simulation.
        forceEqGen = false;  % Force generation of equations.
    end
        
    properties (Dependent)
        thermalVoltage;          % Thermal voltage.
    end
    
    methods
        % constructor
        function obj = circuit(name)
            % CIRCUIT Create an empty circuit.
            % obj = circuit() creates an empty circuit.
            % obj = circuit(cName) creates an empty circuit and names it
            % cName. cName must be a string.
            if nargin > 0
                obj.name = name;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %                    Component handling                           %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function newComp = addComponent(thisCircuit, circComp, varargin)
            % ADDCOMPONENT Add a component to the circuit.
            % newComp = addComponent(comp, varargin) Adds comp to the
            % circuit and returns a handle to the added component in the
            % output argument newComp. The first argument must be an object
            % derived from the abstract class circuitComponent. The
            % following arguments describe the nodes to which the
            % component is to be attached a.k.a. terminals.
            %
            % EXAMPLE: addComponent(resistor(1e3), 'v1', 'gnd');
            % Above is the preferred method to populate a circuit. The
            % nodes are given as strings and are created automatically and
            % added to the circuit if they don't already exist.
            % 
            % Other methods of adding a component are possible. If the only
            % argument to addComponent is a component object, we assume
            % that the component is already connected to the necessary
            % nodes and add those nodes to the circuit. 
            % Alternatively, connection nodes in the argument varargin can
            % be provided as node objects.
            
            % check if the component is already in the circuit
            if any(thisCircuit.components == circComp)
                error('This circuit component is already in the circuit!');
            end
            
            nT = circComp.numTerms;
            
            % provide support for no terminal specification
            if numel(varargin) == 0
                connectionTerms = num2cell(circComp.terminals);
            else
                connectionTerms = varargin;
            end
            
            % check if the correct number of nodes are provided necessary
            % for the connection
            if numel(connectionTerms) ~= nT
                error(['Number of supplied nodes does not match the '...
                       'required terminal number!']);
            end
            
            % loop through all terminals and add them to the circuit in
            % case they are not already in it
            for i = 1:nT
                circComp.setTerminal(i, thisCircuit.addNode(connectionTerms{i}));
            end          
            
            % set the parent circuit of the component to this circuit.
            circComp.setParent(thisCircuit);
            
            % add the circuit component to the list of circuit components
            nComp = thisCircuit.numComponents + 1;
            thisCircuit.numComponents = nComp;
            thisCircuit.components(nComp) = circComp;
            
            % update the number of voltage variables. Only include the
            % internal nodes of the component as new variables.
            thisCircuit.numVoltageVars = thisCircuit.numVoltageVars +...
                                                circComp.numVoltageVars - nT;
            % update the number of current variables
            thisCircuit.numCurrentVars = thisCircuit.numCurrentVars +...
                                                    circComp.numCurrentVars;
            % update the number of voltage sources if necessary
            thisCircuit.numVoltageSources = thisCircuit.numVoltageSources +...
                                                    circComp.numVoltageSources;
            
            % return the component we added to the circuit
            newComp = circComp;
            thisCircuit.isSealed = false;
        end
        
        function comp = findComponent(thisCircuit, aComp)
            % FINDCOMPONENT Find a component in the circuit.
            % Find the reference of aComp in the list of nodes in the
            % circuit and provide its handle. The input aComp can either
            % be a handle or a string.
            %
            % The circuit doesn't have any knowledge of the components of
            % its subcircuits. Hence only top level components can be
            % looked for.
            if ischar(aComp)
                comp = findobj([thisCircuit.components], 'name', aComp);
            elseif isa(aComp, 'circuitComponent')
                comp = thisCircuit.components(thisCircuit.components == aComp);
            else
                error(['You can only search a component with its '...
                     'name or handle!']);
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %                    Node handling                                %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function newNode = addNode(thisCircuit, aNode)
            % ADDNODE Add a node to this circuit.
            % 
            % newNode = addNode(aNode) adds aNode to the circuit. The
            % argument to this method can either be a node object or a
            % string. The added node is returned in the output argument
            % newNode.
            
            % first check if aNode is already in the circuit
            newNode = thisCircuit.findNode(aNode);
            
            % if aNode is not found in the circuit
            if isempty(newNode)
                if isa(aNode, 'node')
                    % a node object is given as argument.
                    if ~isempty(thisCircuit.findNode(aNode.name))
                        error(['There is already a node with the name '...
                                aNode ' in this circuit!']);
                    else
                        newNode = aNode;
                    end
                elseif ischar(aNode)
                    % a string is given as argument. create new node.
                    newNode = node(aNode);
                else
                    error('The argument must be either a node or a string!');
                end
                
                % we do not add nodes without names
                if isempty(newNode.name)
                    error('This node does not have a name!');
                end
                
                % if everything is OK, add the node        
                
                % update the circuit
                nodeNumber = thisCircuit.numNodes + 1;
                thisCircuit.numNodes = nodeNumber;
                thisCircuit.numVoltageVars = thisCircuit.numVoltageVars + 1;
                
                % convey info to the node
                newNode.number = nodeNumber;
                
                % add node to the list of nodes in the circuit
                thisCircuit.nodes(nodeNumber) = newNode;
                newNode.parentCircuit = thisCircuit; 
            end
        end
       
        function n = findNode(thisCircuit, aNode)
            % FINDNODE Find a node in the circuit.
            %   Find the reference of the input node "aNode" in the list of
            %   nodes in the circuit and provide its handle. Search a node
            %   either with its name, number or handle.
            if ischar(aNode)
                n = findobj(thisCircuit.nodes, 'name', aNode);
            elseif isnumeric(aNode)
                n = findobj(thisCircuit.nodes, 'number', aNode);
            elseif isa(aNode, 'node')
                n = thisCircuit.nodes(thisCircuit.nodes == aNode);
            else
                error(['You can only search a node with its '...
                     'name, number or handle!']);
            end
        end
        
        % Ground node handling
        function setGroundNode(thisCircuit, aNode)
            % SETGROUNDNODE set the ground node of the circuit.
            % setGroundNode(aNode) sets aNode (string or node object) as
            % the ground node of the circuit. Note that without a ground
            % node, the KCL equations contain one equation too many and 
            % hence the voltages cannot be determined uniquely.
            n = thisCircuit.findNode(aNode);
            
            if ~isempty(aNode) && isempty(n) 
                error('There is no such node in this Circuit');
            end
            thisCircuit.groundNode = n;
            thisCircuit.Xi = [];
            thisCircuit.Vi = [];
            thisCircuit.ti = 0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %               Sealing of nodes and components                   %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % NOTE: This method should be called from the solver class with
        % the necessary checks before a solution is attempted. I.e.
        % check if the circuit is sealed; if not call this method.
        % Additionally, any operation that causes a reassignment of
        % equation numbers should set thisCircuit.isSealed = false.

        function seal(thisCircuit)
            %SEAL determine the places of variables in the overall
            %equations.
            % SEAL method is to be used after all other operations
            % regarding circuit creation are completed. SEAL assigns
            % numbers to all nodes and components that are required to
            % place their variables to the right locations in the KCL and
            % KVL equations.
            
            if isempty(thisCircuit.groundNode)
                error('Circuit has no ground node!');
            end
            
            % Seal each component
            for i = 1:thisCircuit.numComponents
                cc = thisCircuit.components(i);
                cc.sealComponent;
                % make the circuit aware of the voltage sources
                if cc.numVoltageSources > 0
                    % append component's voltage sources to the list of
                    % voltage sources. If cc is a source device, the
                    % device itself is added to the list.
                    thisCircuit.voltageSourceNodeNums = [...
                                        thisCircuit.voltageSourceNodeNums;
                                        cc.voltEqNums(cc.voltageSourceNodeNums)];
                    thisCircuit.voltageSources = [...
                                                 thisCircuit.voltageSources; 
                                                 cc.voltageSources];
                end
            end
            
            % set the ground node number
            if isempty(thisCircuit.groundNode)
                thisCircuit.groundNodeNumber = 0;
            else
                thisCircuit.groundNodeNumber = ...
                                thisCircuit.groundNode.getTopLevelNumber();
            end
            
            % set number of independent voltage variables
            thisCircuit.numIndepVoltVars = thisCircuit.numVoltageVars -...
                                                ~isempty(thisCircuit.groundNode);
            
            % set total number of variables
            thisCircuit.numVars = thisCircuit.numIndepVoltVars + ...
                                                thisCircuit.numCurrentVars;
                                            
            % set total number of variables for the no currents case
            thisCircuit.numVarsnc = thisCircuit.numIndepVoltVars - ...
                                            thisCircuit.numVoltageSources;
                                        
            thisCircuit.Xi = NaN(thisCircuit.numVars,1);
            thisCircuit.Vi = NaN(thisCircuit.numVarsnc,1);
                   
            % Below we arrage the order of nodes connected to sources.
            % A pair [a,b] in the voltageSourceNodeNums vector (Kx2) means
            % there is a voltage source from node "a" to node "b", where a
            % and b are node numbers.
            %        [....]
            % vsnn = [a, b]
            %        [....]
            % In this case b is a dependent node to a. If a node becomes
            % dependent on two different nodes, we try to reverse the order
            % in the voltageSourceNodeNums entry. If this doesn't work
            % either we exit with an error. This simple check is necessary
            % to be able to handle voltage source chains and this method
            % can detect voltage source loops.
            % FIXME: if there is a voltage source chain
            % A----(V)-----B-----(V)-----C----(V)-----D
            % and if we add the sources in the order
            % A->B, D->C, B->C then our simple algorithm doesn't try to
            % reorder (D->C) and erroneously thinks this is a loop. This
            % situation should be pretty uncommon though.
            vsnn = thisCircuit.voltageSourceNodeNums;
            numVS = size(vsnn,1);
            vsSigns = ones(numVS, 1);
            grndNum = thisCircuit.groundNodeNumber;
            for k = 1:numVS
                % is the dependent node already on the list of dependent
                % nodes?
                if any(vsnn(k, 2) == [vsnn(1:k-1, 2); grndNum])
                    % is the other node dependent to any other node?
                    if any(vsnn(k, 1) == [vsnn(1:k-1, 2); grndNum])
                        error(['This circuit contains a voltage source loop!\n',...
                                'The source between nodes %d -> %d is a ',...
                                'part of it.\n',...
                                'If you have a chain of sources linked to ',...
                                'the ground, add the source directly ',...
                                'connected to the ground before others.'],...
                                vsnn(k,1), vsnn(k,2));
                    end
                    % if not make the other node dependent
                    vsnn(k, [1 2]) = vsnn(k, [2 1]);
                    % reverse the direction of the voltage source
                    vsSigns(k) = -1;
                end
            end
            % update the circuit
            thisCircuit.voltageSourceSigns = vsSigns;
            thisCircuit.voltageSourceNodeNums = vsnn;
            if numVS ~= 0
                thisCircuit.eqNumsnc = setdiff((1:thisCircuit.numVoltageVars)',...
                                                   [grndNum; vsnn(:,2)]);
            else
                thisCircuit.eqNumsnc = setdiff((1:thisCircuit.numVoltageVars)',...
                                                   [grndNum]);
            end
            
            thisCircuit.isSealed = true;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %                    Equation supply functions                    %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
        % Equation generation
        % These methods supply the elements of the circuit equation for a
        % particular input x and time point t. If the variables (x, t) are
        % new to the circuit the equations are generated from scratch.
        % If the circuit already knows the values for x,t they are provided
        % as is. This way device eqations are not recalculated if different
        % parts of the differential equation are needed. Varargin can
        % either be true or false and can be used to force the circuit to
        % regenerate equations.
        
        
        
        function Q = Q(thisCircuit, x, t)
            % Q Get the charge vector Q.
            thisCircuit.statusCheck(x,t);
            Q = thisCircuit.Qi;               
        end
        
        function I = I(thisCircuit, x, t)
            % I Get the current vector I.
            thisCircuit.statusCheck(x,t);
            I = thisCircuit.Ii;               
        end
        
        function J = J(thisCircuit, x, t)
            % J Get the source vector J.
            thisCircuit.statusCheck(x,t);
            J = thisCircuit.Ji;               
        end
        
        function dQ = dQ(thisCircuit, x, t)
            % DQ Get the Jacobian of Q.
            thisCircuit.statusCheck(x,t);
            dQ = thisCircuit.dQi;               
        end
        
        function dI = dI(thisCircuit, x, t)
            % DI Get the Jacobian of I.
            thisCircuit.statusCheck(x,t);
            dI = thisCircuit.dIi;               
        end
        
        function [Q, I, J, dQ, dI] = stamp(thisCircuit, x, t)
            thisCircuit.statusCheck(x,t);
            Q = thisCircuit.Qi;               
            I = thisCircuit.Ii;               
            J = thisCircuit.Ji;               
            dQ = thisCircuit.dQi;               
            dI = thisCircuit.dIi;               
        end
        
       % Functions for the generation of current-free equation set. These
       % equations do not include the voltage source currents explicitly.
        
        function Qnc = Qnc(thisCircuit, v, t)
            % QNC Get the Q vector - no current variables.
            thisCircuit.statusCheckNoCurr(v,t);
            Qnc = thisCircuit.Qinc;               
        end
        
        function Inc = Inc(thisCircuit, v, t)
            % INC Get the I vector - no current variables.
            thisCircuit.statusCheckNoCurr(v,t);
            Inc = thisCircuit.Iinc;               
        end
        
        function Jnc = Jnc(thisCircuit, v, t)
            % JNCGet the J vector - no current variables.
            thisCircuit.statusCheckNoCurr(v,t);
            Jnc = thisCircuit.Jinc;               
        end
        
        function dQnc = dQnc(thisCircuit, v, t)
            % DQNC Get the Jacobian of Q - no current variables.
            thisCircuit.statusCheckNoCurr(v,t);
            dQnc = thisCircuit.dQinc;               
        end
        
        function dInc = dInc(thisCircuit, v, t)
            % DINCGet the Jacobian of I - no current variables.
            thisCircuit.statusCheckNoCurr(v,t);
            dInc = thisCircuit.dIinc;               
        end
        
        function [Qnc, Inc, Jnc, dQnc, dInc] = stampnc(thisCircuit, v, t)
            thisCircuit.statusCheckNoCurr(v,t);
            Qnc = thisCircuit.Qinc;               
            Inc = thisCircuit.Iinc;               
            Jnc = thisCircuit.Jinc;               
            dQnc = thisCircuit.dQinc;               
            dInc = thisCircuit.dIinc;               
        end
        
        function forceEquationGeneration(thisSubcircuit)
            thisSubcircuit.forceEqGen = true;
        end
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %                    Getters / Setters                            %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function vt = get.thermalVoltage(thisCircuit)
            vt = vF.BoltzmannConstant*thisCircuit.Tempi/vF.elementaryCharge;
        end
        
        function Q = getFullQ(thisCircuit, v, t)
            thisCircuit.statusCheckNoCurr(v,t);
            Q = thisCircuit.Qi;
        end
        
        function ncen = getNoCurrEquationNumber(thisCircuit, tlnn)
            % GETNOCURREQUATIONNUMBER variable number of a node - no curr.
            % This method returns the position of a voltage variable
            % corresponding to a node in the entirety of circuit equations.
            
            % here we use the fact that the number of the ground node never
            % appears in the right column of voltageSourceNodeNums!!!
            ind = thisCircuit.voltageSourceNodeNums(:,2);
            gnd = thisCircuit.groundNodeNumber;
            % if the node is a dependent node return 0
            if any(ind == tlnn)
                ncen = -1;
            % if the node is the ground node return 0
            elseif tlnn == gnd
                ncen = 0;
            % else return the equation number
            else
                ncen = tlnn - sum(tlnn > [gnd; ind]);
            end
        end
        
        function en = getEquationNumber(thisCircuit, tlnn)
            % GETEQUATIONNUMBER find the variable for a node in equations.
            % This method returns the position of a voltage variable
            % corresponding to a node in the entirety of circuit equations
            % in the no current variables case.
            gnn = thisCircuit.groundNodeNumber;
            if gnn == 0 || tlnn < gnn
                en = tlnn;
            elseif gnn == tlnn
                en = 0;
            else
                en = tlnn - 1; 
            end
        end
        
        function tlc = getTopLevelCircuit(thisCircuit)
            tlc = thisCircuit;
        end

        function [SV SC] = getSourceVariables(thisCircuit, v, t, Qprev, tprev, varargin)
            % GETSOURCECURRENTS calculate the voltage source currents.
            %
            % Compute source currents for no-current-variable case.
            % if an additional variable qp (the time derivative of the
            % charges) in varargin is given, then the linear charge model
            % is used to calculate the currents. Otherwise the charge
            % function itself is discretized.
            
            % variables for equation system size
            nVV = thisCircuit.numVoltageVars;
            niVV = thisCircuit.numIndepVoltVars;
            nVS = thisCircuit.numVoltageSources;
            eqNums = thisCircuit.eqNumsnc;
            vsnn = thisCircuit.voltageSourceNodeNums;
            vss = thisCircuit.voltageSourceSigns;
            grndNum = thisCircuit.groundNodeNumber;
            
            % if it is the first time point we calculate the currents with
            % the derivative of the voltage vector, vp
            firstTime = (nargin == 6);
            
            % Pad the derivative vector with zeros to make space for the
            % eliminated voltage variables
            vPad = zeros(nVV, 1);
            vPad(eqNums) = v;
            
            lsn = vsnn(:,1); % reference nodes connected to voltage sources
            rsn = vsnn(:,2); % dependent nodes connected to voltage sources
            
            % loop through all voltage sources
            for i = 1:nVS            
                % get the value and time derivative for the voltage source
                V0 = thisCircuit.voltageSources(i).getVoltage(t);
                % set the time value and derivative of the dependent node
                vPad(rsn(i)) = vPad(lsn(i)) - vss(i)*V0;
            end
            SV = vPad(rsn);
            % calculate the left hand side of the KCL equations
           
            vPad(grndNum) = [];
            vPadc = [vPad;zeros(thisCircuit.numCurrentVars,1)];
            I = thisCircuit.I(vPadc, t);
            
            if firstTime == false
                Q = thisCircuit.Q(vPadc, t);
                SC = -1*((Q(1:niVV)-Qprev(1:niVV))/(t-tprev) + I(1:niVV));
            else
                vp = thisCircuit.dQnc(v,t)\(thisCircuit.Jnc(v,t)-...
                                                    thisCircuit.Inc(v,t));
                for i = 1:nVS
                    vpPad = zeros(nVV, 1);
                    vpPad(eqNums) = vp;
                    [~, Vp0] = thisCircuit.voltageSources(i).getVoltage(t);
                    vpPad(rsn(i)) = vpPad(lsn(i)) - vss(i)*Vp0;
                    vpPad(grndNum) = [];
                    dQ = thisCircuit.dQ(vPadc, t);
                end
                SC = -1*(dQ(1:niVV,1:niVV)*vpPad(1:1:niVV) + I(1:niVV));
            end
            
            % expand SC with the ground node
            if grndNum ~= 0
                SC = [SC(1:grndNum-1); 0; SC(grndNum:niVV)];
            end
            % non-zero elements are the source currents 
            SC = SC(rsn);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %                    Subcircuit conversion                        %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function subcirc = subcircuit(thisCircuit, varargin)
            %SUBCIRCUIT convert circuit to a subcircuit
            metaCirc = metaclass(thisCircuit);
            props = metaCirc.PropertyList;
            % create new subcircuit object
            subcirc = subcircuit;
            % copy all properties of the original object to the new one
            for i = 1:numel(props)
                if ~props(i).Dependent
                    propName = props(i).Name;
                    subcirc.(propName) = thisCircuit.(propName);
                end
            end
            
            if nargin > 1
                subcirc.setExternalNodes(varargin{:});
            end
            % NOTE: should we delete the old object?
            %delete(thisCircuit);
        end
    end
    
    methods (Access = protected)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %                    Main equation generation                     %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Equation Handling
        function status = statusCheck(thisCircuit, x, t)
            % STATUSCHECK check if the internal state matches the new
            % state (x,t).
            nV = thisCircuit.numVars;
            if length(x) ~= nV
                error('The variable vector does not have the correct size!');
            end
            
            status = (thisCircuit.ti ~= t || any(thisCircuit.Xi ~= x)...
                                          || thisCircuit.forceEqGen);
            
            if status == true
                thisCircuit.eqGen(x, t);
            end
        end
        
        function status = statusCheckNoCurr(thisCircuit, v, t)
            % STATUSCHECKNOCURR check if the internal state matches the new
            % state (x,t). No currents case.
            nV = thisCircuit.numVarsnc;
            if length(v) ~= nV
                error('The voltage vector does not have the correct size!');
            end
            
            status = (thisCircuit.ti ~= t  || any(thisCircuit.Vi ~= v)...
                                           || thisCircuit.forceEqGen);
            
            if status
                thisCircuit.eqGenNoCurr(v, t);
            end
        end
        
        function eqGen(thisCircuit, x, t)
            % EQGEN generates circuit equation with current variables
            % EQGEN(x, t) generates the KVL and KCL equations of the
            % circuit and stores the charge, current and source vectors in
            % the internal variables Qi, Ii, Ji. Their jacobians are stored
            % in dQi and dIi. The first argument, x, is the vector of
            % node voltages and branch currents (only for necessary
            % branches). If the circuit has a ground node set, the voltage
            % for this node should not be included in x. The voltage for
            % the ground node is appended to x before calling the eqEval
            % method. The second argument is the time. This method is only
            % called by "statusCheck" method if x and t don't match the
            % internal state variables Xi and ti.
            niV = thisCircuit.numIndepVoltVars;

            % Separate the voltage and current variables in x.
            voltageVars = x(1:niV);
            currentVars = x(niV+1:end);
            grndNum = thisCircuit.groundNodeNumber;

            % If there is a ground node, insert a zero in the voltage
            % variables.
            if grndNum == 0
                voltVarsPad = voltageVars;
            else
                voltVarsPad = [voltageVars(1:grndNum-1); 0; voltageVars(grndNum:end)];
            end
            
            xtot = [voltVarsPad; currentVars];
            % call the equation evaluation method
            [Q, I, J, dQ, dI] = thisCircuit.eqEval(xtot, t);
            
            % Delete the rows and columns for the ground node.
            if grndNum ~= 0
                Q(grndNum) = [];
                I(grndNum) = [];
                J(grndNum) = [];
                dQ(grndNum, :) = [];
                dQ(:, grndNum) = [];
                dI(grndNum, :) = [];
                dI(:, grndNum) = [];
            end
            
            % store the state variables
            thisCircuit.Xi = x;
            thisCircuit.ti = t;
            
            % store the state variables for the no current case
            % this is done because getSourceVariables function calls the
            % full equation generation routines in the current free
            % simulations.
            thisCircuit.Vi = xtot(thisCircuit.eqNumsnc);
            
            % store the equations
            thisCircuit.Qi = Q;
            thisCircuit.Ii = I;
            thisCircuit.Ji = J;
            thisCircuit.dQi = dQ;
            thisCircuit.dIi = dI;
            % set force switch to false
            thisCircuit.forceEqGen = false;
        end
        
        function eqGenNoCurr(thisCircuit, v, t)
            % EQGENNOCURR generates circuit equations without current
            % variables.
            % EQGENNOCURR(v, t) generates only the KCL equations of the
            % circuit and additionally eliminates all voltage variables
            % that are dependent to another voltage variable via a voltage
            % source. *THIS ONLY WORKS IF THERE ARE NO INDUCTIVE BRANCHES!*
            % The first argument, v, is the vector of voltages (without the
            % eliminated variables and the ground node voltage).
            %
            % EQGENNOCURR first assigns values to dependent nodes according
            % to corresponding voltage sources and then generates the
            % circuit equations the usual way. In this intermediate step,
            % charges, currents and sources are stored in the internal
            % variables Qi, Ii and Ji respectively. Their jacobians are
            % stored in dQi and dIi. Then the necessary modifications to
            % the equation system is made to eliminate the dependent nodes
            % and the branch currents. The new equations are stored
            % internally in Qinc, Iinc, Jinc, dQinc and dIinc.
            
            %NOTE: this method does not work when there are inductive
            % branches in the circuit. 
            
            nVV = thisCircuit.numVoltageVars;
            nCV = thisCircuit.numCurrentVars;
            nEq = thisCircuit.eqNumsnc;
            nVS = thisCircuit.numVoltageSources;
            vsnn = thisCircuit.voltageSourceNodeNums;
            vsSign = thisCircuit.voltageSourceSigns;
            vsList = thisCircuit.voltageSources;
            
            voltVarsPad = zeros(nVV, 1);
            voltVarsPad(nEq) = v;
            grndNum = thisCircuit.groundNodeNumber;
            eqNumsnc = thisCircuit.eqNumsnc;
            
            lsn = vsnn(:,1); % reference nodes connected to voltage sources
            rsn = vsnn(:,2); % dependent nodes connected to voltage sources
            
            for i=1:nVS
                n1 = lsn(i); 
                n2 = rsn(i);
                voltVarsPad(n2) = voltVarsPad(n1) - vsSign(i)*(vsList(i).getVoltage(t));
            end
            
            % generate standard equations with current variables
            [Q, I, J, dQ, dI] = thisCircuit.eqEval([voltVarsPad; zeros(nCV, 1)], t);

            
            % store the full equations
            thisCircuit.Xi = [voltVarsPad; zeros(nCV,1)];
            thisCircuit.Qi = Q;
            thisCircuit.Ii = I;
            thisCircuit.Ji = J;
            thisCircuit.dQi = dQ;
            thisCircuit.dIi = dI;
            
            thisCircuit.Xi(grndNum) = [];
            thisCircuit.Qi(grndNum) = [];
            thisCircuit.Ii(grndNum) = [];
            thisCircuit.Ji(grndNum) = [];
            thisCircuit.dQi(grndNum,:) = [];
            thisCircuit.dQi(:,grndNum) = [];
            thisCircuit.dIi(grndNum,:) = [];
            thisCircuit.dIi(:,grndNum) = [];
            
            % truncate the equation system
            Q = Q(1:nVV);
            I = I(1:nVV);
            J = J(1:nVV);
            dQ = dQ(1:nVV,1:nVV);
            dI = dI(1:nVV,1:nVV);            
            
            % add the dependent rows/columns to reference rows/columns
            Q(lsn) = Q(lsn) + Q(rsn);
            I(lsn) = I(lsn) + I(rsn);
            dQ(lsn,:) = dQ(lsn,:) + dQ(rsn,:);
            dQ(:,lsn) = dQ(:,lsn) + dQ(:,rsn);
            dI(lsn,:) = dI(lsn,:) + dI(rsn,:);
            dI(:,lsn) = dI(:,lsn) + dI(:,rsn);

            % store the internal state
            thisCircuit.Vi = v;
            thisCircuit.ti = t;
            
            % Delete the rows and columns for the ground node and dependent
            % variables. No problem if grndNum is repeated in rsn.
            % store reduced equations
            thisCircuit.Qinc = Q(eqNumsnc);
            thisCircuit.Iinc = I(eqNumsnc);
            thisCircuit.Jinc = J(eqNumsnc);
            thisCircuit.dQinc = dQ(eqNumsnc, eqNumsnc);
            thisCircuit.dIinc = dI(eqNumsnc, eqNumsnc);

            % set force switch to false
            thisCircuit.forceEqGen = false;
        end

        function [Q I J dQ dI] = eqEval(thisCircuit, x, t)
            % EQEVAL Evaluate the individual parts of the circuit equation
            % 
            % EQEVAL(x, t) evaluates the parts of the KCL and KVL equations
            % by iterating over all components in the circuit. The first
            % argument, x, is a Nx1 column vector, where N is the number of
            % nodes plus the number of current variables (these usually
            % stem from voltage sources and inductors). This method does
            % not handle any ground node specific function.
            %
            % This method should be overloaded by any compact device model
            % that is set up as a subcircuit. For an example see the
            % MOSFET models included.
            
            nVV = thisCircuit.numVoltageVars;
            nEq = thisCircuit.numVoltageVars + thisCircuit.numCurrentVars;
            
            % Allocate space for equation vectors and matrices.
            Q = zeros(nEq, 1);
            I = zeros(nEq, 1);
            J = zeros(nEq, 1);
            dQ = spalloc(nEq, nEq, nEq*5);
            dI = spalloc(nEq, nEq, nEq*5);
            
            % Separate the voltage and current variables in x.
            voltageVars = x(1:nVV);
            currentVars = x(nVV+1:end);
            
            for i=1:thisCircuit.numComponents
                circComp = thisCircuit.components(i);
                vNums = circComp.voltEqNums;
                cNums = circComp.currEqNums;
                
                v = voltageVars(vNums);
                c = currentVars(cNums);
                totVars = [v; c];
                
                [sQ, sI, sJ, sdQ, sdI] = circComp.stamp(totVars, t);
                  
                eqNums = [vNums (nVV+cNums)];
            
                % Put the entries of the component to their respective
                % places in the overall circuit equations.
                Q = Q + sparse(eqNums,1,sQ,nEq,1);
                I = I + sparse(eqNums,1,sI,nEq,1);
                J = J + sparse(eqNums,1,sJ,nEq,1);
                
                [idQ, jdQ, sdQ] = find(sdQ);
                [idI, jdI, sdI] = find(sdI);
                
                dQ = dQ + sparse(eqNums(idQ), eqNums(jdQ), sdQ, nEq, nEq);
                dI = dI + sparse(eqNums(idI), eqNums(jdI), sdI, nEq, nEq);
                
                if isa(circComp, 'neuronHH') == true
                    if circComp.numExReceptors > 0
                        for k=1:1:circComp.numExReceptors
                            s = circComp.getCrossRates(circComp.exSourceNeuron{k}, totVars, k, 'exReceptor');
                            dI(cNums(3+k),circComp.exSourceNeuron{k}.voltEqNums) = dI(cNums(3+k),circComp.exSourceNeuron{k}.voltEqNums) + s;
                        end
                    end
                    
                    if circComp.numInReceptors > 0
                        for k=1:1:circComp.numInReceptors
                            s = circComp.getCrossRates(circComp.inSourceNeuron{k}, totVars, k, 'inReceptor');
                            dI(cNums(3+circComp.numExReceptors+k),circComp.inSourceNeuron{k}.voltEqNums) = dI(cNums(3+circComp.numExReceptors+k),circComp.inSourceNeuron{k}.voltEqNums) + s;
                        end
                    end
                end
                
            end
        end        
        
        
    end
    
    methods (Static)
        function tlnn = getTopLevelNodeNumber(aNodeNumber)
            % GETTOPLEVELNODENUMBER find the location of a variable.
            % This method can be used to find the variable number in the
            % overall circuit equations that corresponds to a node. 
            % NOTE: this method does not care about ground node or
            % eliminated nodes in the no-current case. The return number of
            % this method must be processed further according to the type
            % of equations to get the actual equation number.
            tlnn = aNodeNumber;
        end
    end
end

