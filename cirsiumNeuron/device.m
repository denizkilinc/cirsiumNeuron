classdef device < circuitComponent
    % DEVICE Abstract class to serve as a base for devices.
    % This class defines the methods a device must implement, but also
    % implements a lot of the functionality of a device itself. It leaves
    % only the basic device properties for individual device classes such
    % as resistor to implement.
    % The equation of the whole circuit is assumed to be formulated as:
    %                   Q'(v) + I(v) = J
    
    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------

    properties (Abstract, Constant)
        numTerms;
        numCurrentVars;
    end
    
    methods (Abstract)
        sQ = Q(thisCircuitComponent, x, t); % Q vector of the device.
        sI = I(thisCircuitComponent, x, t); % I vector of the device.
        sJ = J(thisCircuitComponent, x, t); % J vector of the device.
        sdQ = dQ(thisCircuitComponent, x, t); % Jacobian of the Q vector.
        sdI = dI(thisCircuitComponent, x, t); % Jacobian of the I vector.
    end
    
    methods 
        function obj = device(name)
            % DEVICE constructor for a device.
            %   All arguments are optional. A real device can be created
            %   with the constructor of the child class.
            if nargin > 0
                obj.name = name;
            end
            obj.numVoltageVars = obj.numTerms;
            obj.numVars = obj.numTerms + obj.numCurrentVars;
        end       
        
        function [Q, I, J, dQ, dI] = stamp(thisDevice, x, t)
            Q = thisDevice.Q(x,t);
            I = thisDevice.I(x,t);
            J = thisDevice.J(x,t);
            dQ = thisDevice.dQ(x,t);
            dI = thisDevice.dI(x,t);
        end
        
        function setEquationNumbers(thisDevice)
            currInd = thisDevice.currentVarInd;
            
            vNums = thisDevice.termNodeNumbers;
            cNums = currInd:currInd+thisDevice.numCurrentVars-1;
            
            thisDevice.voltEqNums = vNums;
            thisDevice.currEqNums = cNums;
        end
        
        function tlc = getTopLevelCircuit(thisDevice)
            if isempty(thisDevice.parentCircuit)
                error('This device does not have a parent circuit')
            else
                tlc = thisDevice.parentCircuit.getTopLevelCircuit();
            end
        end

    end
    
    methods (Access = protected)
        function newDev = copyElement(thisDevice)
            newDev = copyElement@circuitComponent(thisDevice);
            newDev.numVoltageVars = thisDevice.numVoltageVars;
        end
    end

  
    
end

