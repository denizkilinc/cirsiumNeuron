% Simulation setup for the results presented in Section 6.3.
% Analysis of fully connected EI network with heterogeneous neurons via steady-state 
% semi-analytical method.

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

membraneArea = 1000e-12;
currentDensity = 25e-2;
T_vIn = 400e-3;
Ihigh = currentDensity*membraneArea;

neuron1 = neuronMC(2*membraneArea, 0, 4, 4, 'HH', 'Neuron1');
neuron2 = neuronMC(8*membraneArea, 0, 4, 4, 'HH', 'Neuron2');
neuron3 = neuronMC(4*membraneArea, 0, 4, 4, 'HH', 'Neuron3');
neuron4 = neuronMC(1*membraneArea, 0, 4, 4, 'HH', 'Neuron4');
neuron5 = neuronMC(5*membraneArea, 0, 4, 4, 'HH', 'Neuron5');

f_vIn = @(t)(cosTapRect(t-0e-3, T_vIn, 1/1000, 1, 1));
iIn1 = currentSource(2*Ihigh, f_vIn, 'Iin1');
iIn2 = currentSource(8*Ihigh, f_vIn, 'Iin2');
iIn3 = currentSource(4*Ihigh, f_vIn, 'Iin3');
iIn4 = currentSource(1*Ihigh, f_vIn, 'Iin4');
iIn5 = currentSource(5*Ihigh, f_vIn, 'Iin5');

neuron5.synapse(neuron1,'exSlowReceptorMC');
neuron1.synapse(neuron2,'exSlowReceptorMC');
neuron2.synapse(neuron3,'exSlowReceptorMC');
neuron3.synapse(neuron4,'exSlowReceptorMC');
neuron4.synapse(neuron5,'exSlowReceptorMC');

neuron2.synapse(neuron1,'exSlowReceptorMC');
neuron3.synapse(neuron2,'exSlowReceptorMC');
neuron4.synapse(neuron3,'exSlowReceptorMC');
neuron5.synapse(neuron4,'exSlowReceptorMC');
neuron1.synapse(neuron5,'exSlowReceptorMC');

neuron3.synapse(neuron1,'exSlowReceptorMC');
neuron4.synapse(neuron2,'exSlowReceptorMC');
neuron5.synapse(neuron3,'exSlowReceptorMC');
neuron1.synapse(neuron4,'exSlowReceptorMC');
neuron2.synapse(neuron5,'exSlowReceptorMC');

neuron4.synapse(neuron1,'exSlowReceptorMC');
neuron5.synapse(neuron2,'exSlowReceptorMC');
neuron1.synapse(neuron3,'exSlowReceptorMC');
neuron2.synapse(neuron4,'exSlowReceptorMC');
neuron3.synapse(neuron5,'exSlowReceptorMC');


neuron5.synapse(neuron1,'inReceptorMC');
neuron1.synapse(neuron2,'inReceptorMC');
neuron2.synapse(neuron3,'inReceptorMC');
neuron3.synapse(neuron4,'inReceptorMC');
neuron4.synapse(neuron5,'inReceptorMC');

neuron2.synapse(neuron1,'inReceptorMC');
neuron3.synapse(neuron2,'inReceptorMC');
neuron4.synapse(neuron3,'inReceptorMC');
neuron5.synapse(neuron4,'inReceptorMC');
neuron1.synapse(neuron5,'inReceptorMC');

neuron3.synapse(neuron1,'inReceptorMC');
neuron4.synapse(neuron2,'inReceptorMC');
neuron5.synapse(neuron3,'inReceptorMC');
neuron1.synapse(neuron4,'inReceptorMC');
neuron2.synapse(neuron5,'inReceptorMC');

neuron4.synapse(neuron1,'inReceptorMC');
neuron5.synapse(neuron2,'inReceptorMC');
neuron1.synapse(neuron3,'inReceptorMC');
neuron2.synapse(neuron4,'inReceptorMC');
neuron3.synapse(neuron5,'inReceptorMC');

ckt = circuitMC('Single Neuron Test');

ckt.addComponent(neuron1, 'Node1', 'gnd');
ckt.addComponent(neuron2, 'Node2', 'gnd');
ckt.addComponent(neuron3, 'Node3', 'gnd');
ckt.addComponent(neuron4, 'Node4', 'gnd');
ckt.addComponent(neuron5, 'Node5', 'gnd');
ckt.addComponent(iIn1, 'gnd', 'Node1');
ckt.addComponent(iIn2, 'gnd', 'Node2');
ckt.addComponent(iIn3, 'gnd', 'Node3');
ckt.addComponent(iIn4, 'gnd', 'Node4');
ckt.addComponent(iIn5, 'gnd', 'Node5');
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
