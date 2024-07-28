% Simulation setup for the results presented in Section 6.2.
% Analysis of single spiking neuron via steady-state semi-analytical method.

%-------------------------------------------------------------------------------------
% Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
% All Rights Reserved
%-------------------------------------------------------------------------------------

clear classes;
clear all
close all
clc

addpath('../cirsiumNeuron');
addpath('../cirsiumNeuron/steady_state_method');
addpath('../cirsiumNeuron/transient_method');

membraneArea = 5*1000e-12;
currentDensity = 25e-2;
T_vIn = 1000e-3;
Ihigh = currentDensity*membraneArea;

neuron1 = neuronMC(membraneArea, 0, 0, 0, 'HH', 'Neuron1');

f_vIn = @(t)(cosTapRect(t-0e-3, T_vIn, 1/1000, 1, 1));
iIn1 = currentSource(Ihigh, f_vIn, 'Iin1');

ckt = circuitMC('Single Neuron Test');

ckt.addComponent(neuron1, 'Node1', 'gnd');
ckt.addComponent(iIn1, 'gnd', 'Node1');
ckt.setGroundNode('gnd');
ckt.seal();

% solver settings
solverType = 'SEQ-TRPZ';
linSolverType = 'GMRES';
currentVariables = 'on';
showSolution = 'off';
t0 = 0;
tend = 1*T_vIn;
timeStep = 1e-5;

% initial conditions
nTr = ckt.numMCs;
switch currentVariables
    case 'on'
        nV = ckt.numVars;
        nVV = ckt.numIndepVoltVars;
    case 'off'
        nV = ckt.numVarsnc;
        nVV = nV;
end

y0 = -0.065*ones(nV,1); %initial state vector for membrane voltages
s0 = [];                %initial state vector for ion channels

inhibitory_channels = 20; %number inhibitory receptor channels per synapse
excitatory_channels = 20; %number excitatory receptor channels per synapse

for i=1:1:nTr
    if isa(ckt.MCs{i},'inReceptorMC') == true
        s0 = [s0; inhibitory_channels*ckt.MCs{i}.stateVector];
    elseif isa(ckt.MCs{i},'exSlowReceptorMC') == true
        s0 = [s0; excitatory_channels*ckt.MCs{i}.stateVector];
    else
        s0 = [s0; ckt.MCs{i}.stateVector*(membraneArea/1e-12)];
    end
end

%% noiseless transient simulation
solver = transientSolverMC(ckt, 'solverType', solverType,...
                            'linSolverType', linSolverType,...
                            'showSolution', showSolution,... 
                            'currentVariables', currentVariables,...
                            'y0', y0,...
                            'yp0', [],...
                            't0', t0,...
                            's0', s0,...
                            'sp0', [],...
                            'seq0', s0,...
                            'timeStep', timeStep,...
                            'tend', tend,...
                            'breakPoints', [],...
                            'reltol', 1e-3,...
                            'abstol_v', 1e-6,...
                            'abstol_c', 1e-9,...
                            'abstol_q', 1e-21,...
                            'chargeScaleFactor', 1e4/T_vIn,...
                            'lmax', 1e12);

tic
solver.solve();
duration = toc;

solver.displaySolution;

%% estimate steady state solution with transient sim
Tf_est = 0.0108; %this value should be obtained from noiseless transient simulation 
t0 = tend - 2*Tf_est;
hbal = harmonicBalanceSolverMC(ckt, 'nFreq', 100,...
                                  'Tf_est', Tf_est,...
                                  't0', t0,...
                                  'yIdx', 1,...
                                  'currentVariables', currentVariables,...
                                  'reltol', 1e-3,...
                                  'abstol_v', 1e-6,...
                                  'abstol_c', 1e-9,...
                                  'linSolverType', 'GMRES',...
                                  'gmrestart', 250,...
                                  'gpufft', false);
hbal.estimateSS(solver);

tic
hbal.solve();
toc

Tf = hbal.Tf; %oscillation period of the neuronal circuit

hbal.displaySolution;

%% LPTV
sys = lptvMC;
sys.lptvGen(hbal);

c_ss = sys.computeJS(1); %asymptotic timing jitter variance slope
