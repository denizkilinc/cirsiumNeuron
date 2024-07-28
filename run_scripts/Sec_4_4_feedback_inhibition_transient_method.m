% Simulation setup for the results presented in Section 4.4.
% Analysis of feedback inhibition via transient semi-analytical method.

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
T_vIn = 200e-3;
Ihigh = currentDensity*membraneArea;

neuron1 = neuronMC(membraneArea, 0, 0, 1, 'MC', 'Neuron1');
neuron2 = neuronMC(membraneArea, 0, 1, 0, 'MC', 'Neuron2');

neuron2.synapse(neuron1,'inReceptorMC');
neuron1.synapse(neuron2,'exSlowReceptorMC');

f_vIn = @(t)(cosTapRect(t-110e-3, T_vIn, 1/100, 2/200, 1));
iIn = currentSource(Ihigh, f_vIn, 'Iin');
ckt = circuitMC('Single Neuron Test');
ckt.addComponent(neuron1, 'Node1', 'gnd');
ckt.addComponent(neuron2, 'Node2', 'gnd');
ckt.addComponent(iIn, 'gnd', 'Node1');
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

y0 = [-0.065;-0.065]; %initial state vector for membrane voltages
s0=[];                %initial state vector for ion channels

inhibitory_channels = 20000; %number inhibitory receptor channels per synapse
excitatory_channels = 30000; %number excitatory receptor channels per synapse

for i=1:1:nTr
    if isa(ckt.MCs{i},'inReceptorMC') == true
        s0 = [s0; inhibitory_channels*ckt.MCs{i}.stateVector];
    elseif isa(ckt.MCs{i},'exSlowReceptorMC') == true
        s0 = [s0; excitatory_channels*ckt.MCs{i}.stateVector];
    else
        s0 = [s0; ckt.MCs{i}.stateVector*(membraneArea/1e-12)];
    end
end

%% noisy transient simulation
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

%% variance calculation with transient method
solverType = 'TRPZ-RK';
stateSpaceType = 'standard';
timeStep = 2.0e-6;
tref = tend;
Kind = [1,2];

solverTD = nonMonteCarloSolverMC(ckt, 'solverType', solverType,...
                            'linSolverType', linSolverType,...
                            'showSolution', showSolution,... 
                            'currentVariables', currentVariables,...
                            'y0', y0,...
                            'yp0', [],...
                            't0', t0,...
                            's0',s0,...
                            'sp0', [],...
                            'seq0', s0,...
                            'timeStep', timeStep,...
                            'tend', tend,...
                            'tref', tref,...
                            'Kind', Kind,...
                            'breakPoints', [],...
                            'reltol', 1e-3,...
                            'abstol_v', 1e-6,...
                            'abstol_c', 1e-9,...
                            'abstol_q', 1e-21,...
                            'chargeScaleFactor', 1e2/T_vIn,...
                            'lmax', 1e12);
tic
solverTD.solve();
toc

solverTD.displaySolution;

%%
%The simulation duration is 200msec. For the first 100ms, we wait for the ion
%channels in the neurons to reach steady state. Then, a transient current
%stimulus with a duration of 2msec is injected at time t=110msec.
figure
semilogy(solverTD.T(50000:end),solverTD.K(1,50000:end)); grid on; xlim([0.1,0.2]);
title('Var[V_1(t)]');

