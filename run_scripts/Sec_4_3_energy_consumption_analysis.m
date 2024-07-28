% Simulation setup for the results presented in Section 4.3.
% Energy consumption analysis of feedback inhibition.

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
T_vIn = 1000e-3;
Ihigh = currentDensity*membraneArea;

neuron1 = neuronMC(membraneArea, 0, 0, 1, 'MC', 'Neuron1');
neuron2 = neuronMC(membraneArea, 0, 1, 0, 'MC', 'Neuron2');

neuron2.synapse(neuron1,'inReceptorMC');
neuron1.synapse(neuron2,'exSlowReceptorMC');

f_vIn = @(t)(cosTapRect(t-0e-3, T_vIn, 1/1000, 1, 1));
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

inhibitory_channels = 10000; %number inhibitory receptor channels per synapse
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

%% power & energy consumption analysis
neuron = solver.circuit.components(1);

E_K = neuron.E_K;
E_Na = neuron.E_Na;
E_NMDA = neuron.E_NMDA;
E_Cl = neuron.E_Cl;
E_L = neuron.E_L;

g_K = neuron.g_K;
g_Na = neuron.g_Na;
g_ExSlowReceptor = neuron.g_ExSlowReceptor;
g_InReceptor = neuron.g_InReceptor;
g_L = neuron.g_L;

Y_1 = solver.Y(1,:);
Y_2 = solver.Y(2,:);
T = solver.T;

%neuron1
N_K_1 = solver.S(5,:);
N_Na_1 = solver.S(13,:);
N_In_1_1 = solver.S(17,:);
N_In_1_2 = solver.S(18,:);

G_K_1 = g_K*N_K_1;
G_Na_1 = g_Na*N_Na_1;
G_In_1_1 = g_InReceptor*N_In_1_1;
G_In_1_2 = g_InReceptor*N_In_1_2;
G_L_1 = g_L*neuron.Area;

power_K_1 = G_K_1.*(Y_1-E_K).^2;
power_Na_1 = G_Na_1.*(Y_1-E_Na).^2;
power_In_1_1 = G_In_1_1.*(Y_1-E_Cl).^2;
power_In_1_2 = G_In_1_2.*(Y_1-E_Cl).^2;
power_L_1 = G_L_1.*(Y_1-E_L).^2;

power_total_1 = power_K_1 + power_Na_1 + power_In_1_1 + power_In_1_2 + power_L_1;

%neuron2
N_K_2 = solver.S(23,:);
N_Na_2 = solver.S(31,:);
N_ExSlow_2 = solver.S(36,:);

G_K_2 = g_K*N_K_2;
G_Na_2 = g_Na*N_Na_2;
G_ExSlow_2 = g_InReceptor*N_ExSlow_2;
G_L_2 = g_L*neuron.Area;

power_K_2 = G_K_2.*(Y_2-E_K).^2;
power_Na_2 = G_Na_2.*(Y_2-E_Na).^2;
power_ExSlow_2 = G_ExSlow_2.*(Y_2-E_NMDA).^2;
power_L_2 = G_L_2.*(Y_2-E_L).^2;

power_total_2 = power_K_2 + power_Na_2 + power_ExSlow_2 + power_L_2;

power_source = solver.circuit.components(3).A.*Y_1;
power_total = power_total_1 + power_total_2;

%%
figure
plot(T,power_total(1,:));

energy_total = cumtrapz(T,power_total(1,:));

figure
plot(T,energy_total);

%avg_energy = (energy_total(floor(0.1898/(T(2)-T(1)))) - energy_total(floor(0.08185/(T(2)-T(1)))))/7;

%%
figure
plot(T,power_source);

energy_total_source = cumtrapz(T,power_source);

figure
plot(T,energy_total_source);

%avg_energy_source = (energy_total_source(floor(0.1898/(T(2)-T(1)))) - energy_total_source(floor(0.08185/(T(2)-T(1)))))/7;


