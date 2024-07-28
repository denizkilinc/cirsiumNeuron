classdef electrical < handle
    % ELECTRICAL Top level class for circuits.
    
    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------
    
    properties (SetAccess = protected)
        name;                   % Name of the device.
        numVars = 0;            % Number of total variables
        groundNode = node.empty;% Ground node of the circuit.
        groundNodeNumber = 0;
        components = resistor.empty;
        numComponents = 0;      % Number of circuit components.
        numVoltageVars = 0;     % Number of voltage variables.
        numVoltageSources = 0;  % Number of voltage sources.
        voltageSourceNodeNums = [];  % K x 2 array holding the V.S. node numbers.
        voltageSources = voltageSource.empty; % List of independent voltage sources.
        numMCs = 0;             % Number of MCs.
        MCs;
        circCompMCs;
        numCircCompMCs = 0;
        numStates = 0;
        numTransitions = 0;
        isSealed = false;       
        timeDependent = false;
    end
    
    methods (Abstract)
        sQ = Q(thisCircuitComponent, x, t); % Q vector of the device.
        sI = I(thisCircuitComponent, x, t); % I vector of the device.
        sJ = J(thisCircuitComponent, x, t); % J vector of the device.
        sdQ = dQ(thisCircuitComponent, x, t); % Jacobian of the Q vector.
        sdI = dI(thisCircuitComponent, x, t); % Jacobian of the I vector.
        [Q, I, J, dQ, dI] = stamp(thisDevice, x, t)
    end
end

