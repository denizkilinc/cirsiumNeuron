classdef neuronMC < circuitComponentMC
    % NEURONMC Container and controller class for a neuron with ion channels.

    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------
    
    properties (Constant)
        E_K  = -77*10^-3;
        E_Na = 50*10^-3;
        E_L  = -54.387*10^-3;
        E_Cl = -70*10^-3;
        E_AMPA = 0;
        E_NMDA = -20*10^-3;
        g_K  = 20*10^-12;           % single K channel conductance
        g_Na = 20*10^-12;           % single Na channel conductance
        g_L  = 3;                   % total leakage conductance per m^2
        g_InReceptor = 20*10^-12;
        g_ExReceptor = 20*10^-12;
        g_ExSlowReceptor = 20*10^-12;
        c_m  = 10^-2;
        
    end
    
    properties (SetAccess = protected)
        numKChannels = 1;    % number of K Channels
        numNaChannels = 1;   % number of Na Channels
        numInReceptors = 0;  % number of inhibitory receptors
        numExReceptors = 0;  % number of excitatory receptors
        numExSlowReceptors = 0;  % number of slow excitatory receptors
        type = 'MC';
        Area = 0;
        totalKConductance = 0;
        totalNaConductance = 0;
        totalInReceptorConductance = 0;
        totalExReceptorConductance = 0;
        totalExSlowReceptorConductance = 0;
        
    end
    
    methods
        function thisNeuron = neuronMC(Area, numExInputs, numExSlowInputs, numInInputs, type, varargin)
            thisNeuron = thisNeuron@circuitComponentMC(varargin{:});
            thisNeuron.addNode('A');
            thisNeuron.addNode('B');
            thisNeuron.setExternalNodes('A','B');
            thisNeuron.Area = Area;
            thisNeuron.numExReceptors = numExInputs;
            thisNeuron.numExSlowReceptors = numExSlowInputs;
            thisNeuron.numInReceptors = numInInputs;
            thisNeuron.type = type;
            thisNeuron.initMCs();
            
        end
        
        function [rateVector, dRateVector, dRateVectorState, dQst, dIst] = MCRates(thisNeuron, aChannel)
            
            % thisNeuronCh.Xi keeps the node voltages of thisNeuronCh.
            % Xi: Array of voltage and current variables in the circuit.
            % Vi: Array of voltage variables to be used in current-free equations.
            
            
            if isa(aChannel,'KChannelMC') == true
                
                x(1) = thisNeuron.Xi(1);
                x(2) = thisNeuron.Xi(2);
                
                alpha_n = (0.01*((x(1)-x(2))*10^3+55))*10^3/(1-exp(-0.1*((x(1)-x(2))*10^3+55)));
                beta_n  = 0.125*10^3*exp(-0.0125*((x(1)-x(2))*10^3+65));
                
                rateVector=[4*alpha_n; beta_n; 3*alpha_n; 2*beta_n; 2*alpha_n; 3*beta_n; alpha_n; 4*beta_n];
                
                dAlpha_n = ((10^4*(1-exp(-0.1*((x(1)-x(2))*10^3+55)))-10^3*((x(1)-x(2))*10^3+55)*exp(-0.1*((x(1)-x(2))*10^3+55)))*((1-exp(-0.1*((x(1)-x(2))*10^3+55)))^-2));
                dBeta_n = (-1562.5*exp(-0.0125*((x(1)-x(2))*10^3+65)));
                
                dRateVector(:,1)=[4*dAlpha_n; dBeta_n; 3*dAlpha_n; 2*dBeta_n; 2*dAlpha_n; 3*dBeta_n; dAlpha_n; 4*dBeta_n];
                dRateVector(:,2)=-dRateVector(:,1);
                
                dRateVectorState=zeros(aChannel.numTransitions,thisNeuron.numStates);
                
            elseif isa(aChannel,'NaChannelMC') == true
                
                x(1) = thisNeuron.Xi(1);
                x(2) = thisNeuron.Xi(2);
                
                alpha_m = (0.1*((x(1)-x(2))*10^3+40))*10^3/(1-exp(-0.1*((x(1)-x(2))*10^3+40)));
                beta_m  = 4*10^3*exp(-0.0556*((x(1)-x(2))*10^3+65));
                alpha_h = 0.07*10^3*exp(-0.05*((x(1)-x(2))*10^3+65));
                beta_h  = 10^3/(1+exp(-0.1*((x(1)-x(2))*10^3+35)));
                
                rateVector=[3*alpha_m; beta_m; 2*alpha_m; 2*beta_m; alpha_m; 3*beta_m; 3*alpha_m; beta_m; 2*alpha_m; 2*beta_m; alpha_m; 3*beta_m; alpha_h; beta_h; alpha_h; beta_h; alpha_h; beta_h; alpha_h; beta_h];
                
                dAlpha_m = ((10^5*(1-exp(-0.1*((x(1)-x(2))*10^3+40)))-10^4*((x(1)-x(2))*10^3+40)*exp(-0.1*((x(1)-x(2))*10^3+40)))*((1-exp(-0.1*((x(1)-x(2))*10^3+40)))^-2));
                dBeta_m = (-0.2224*10^6*exp(-0.0556*((x(1)-x(2))*10^3+65)));
                
                dAlpha_h = (-35*10^2*exp(-0.05*((x(1)-x(2))*10^3+65)));
                dBeta_h = ((1+exp(-0.1*((x(1)-x(2))*10^3+35)))^-2*10^5*exp(-0.1*((x(1)-x(2))*10^3+35)));
                
                dRateVector=[3*dAlpha_m; dBeta_m; 2*dAlpha_m; 2*dBeta_m; dAlpha_m; 3*dBeta_m; 3*dAlpha_m; dBeta_m; 2*dAlpha_m; 2*dBeta_m; dAlpha_m; 3*dBeta_m; dAlpha_h; dBeta_h; dAlpha_h; dBeta_h; dAlpha_h; dBeta_h; dAlpha_h; dBeta_h];
                dRateVector(:,2)=-dRateVector(:,1);
                
                dRateVectorState = zeros(aChannel.numTransitions,thisNeuron.numStates);
                
            elseif isa(aChannel,'inReceptorMC') == true
                
                x(1) = aChannel.sourceNeuron.Xi(1);
                x(2) = aChannel.sourceNeuron.Xi(2);
                
                Rb1=20*10^6;
                Rb2=10*10^6;
                Ru1=4.6*10^3;
                Ru2=9.2*10^3;
                Ro1=3.3*10^3;
                Ro2=10.6*10^3;
                Rc1=9.8*10^3;
                Rc2=410;
                
                Vpre=x(1)-x(2);
                L=(2.84*10^-3)/(1+exp(-(Vpre*10^3-2)/5));
                
                rateVector=[Rb1*L; Ru1; Rb2*L; Ru2; Ro1; Rc1; Ro2; Rc2];
                dRateVector=zeros(aChannel.numTransitions,2);
                dRateVectorState = zeros(aChannel.numTransitions,thisNeuron.numStates);
                
            elseif isa(aChannel,'exReceptorMC') == true
                
                x(1) = aChannel.sourceNeuron.Xi(1);
                x(2) = aChannel.sourceNeuron.Xi(2);
                
                Rb=13*10^6;
                Ru1=5.9;
                Ru2=8.6*10^4;
                Rd=900;
                Rr=64;
                Ro=2.7*10^3;
                Rc=200;
                
                Vpre=x(1)-x(2);
                L=(2.84*10^-3)/(1+exp(-(Vpre*10^3-2)/5));
                
                rateVector=[Rb*L; Ru1; Rb*L; Ru2; Ro; Rc; Rd; Rr; Rd; Rr];
                dRateVector=zeros(aChannel.numTransitions,2);
                dRateVectorState = zeros(aChannel.numTransitions,thisNeuron.numStates);
                
           elseif isa(aChannel,'exSlowReceptorMC') == true
                
                x(1) = aChannel.sourceNeuron.Xi(1);
                x(2) = aChannel.sourceNeuron.Xi(2);
                
                Rb=1*10^6;
                Ru1=30.5;
                Ru2=16*10^4;
                Rd=840;
                Rr=34;
                Ro=930;
                Rc=369;
                
                Vpre=x(1)-x(2);
                L=(2.84*10^-3)/(1+exp(-(Vpre*10^3-2)/5));
                                
                rateVector=[Rb*L; Ru1; Rb*L; Ru2; Rd; Rr; Ro; Rc];
                dRateVector=zeros(aChannel.numTransitions,2);
                dRateVectorState = zeros(aChannel.numTransitions,thisNeuron.numStates);
                
            else
                
                error('This channel is not a valid channel object!');
                
            end
            
            x(1) = thisNeuron.Xi(1);
            x(2) = thisNeuron.Xi(2);
                
            dQst = zeros(2,thisNeuron.numStates);
            dIst = zeros(2,thisNeuron.numStates);
            
            dIst(1,5) = (x(1)-x(2)-thisNeuron.E_K)*thisNeuron.g_K;
            dIst(2,5) = -dIst(1,5);
            
            dIst(1,13) = (x(1)-x(2)-thisNeuron.E_Na)*thisNeuron.g_Na;
            dIst(2,13) = -dIst(1,13);
            
            numIn = thisNeuron.numInReceptors;
            numEx = thisNeuron.numExReceptors;
            numExSlow = thisNeuron.numExSlowReceptors;
            
            for i =1:1:numIn
                dIst(1,13+i*5-1) = (x(1)-x(2)-thisNeuron.E_Cl)*thisNeuron.g_InReceptor;
                dIst(1,13+i*5) = (x(1)-x(2)-thisNeuron.E_Cl)*thisNeuron.g_InReceptor;
                dIst(2,13+i*5-1:13+i*5) = -dIst(1,13+i*5-1:13+i*5);
            end
            
            for i =1:1:numEx
                dIst(1,13+numIn*5+i*6) = (x(1)-x(2)-thisNeuron.E_AMPA)*thisNeuron.g_ExReceptor;
                dIst(2,13+numIn*5+i*6) = -dIst(1,13+numIn*5+i*6);
            end
            
            B = 1;%1/(1+exp(-0.062*1e3*(x(1)-x(2)))*(1.5/3.57));
            for i =1:1:numExSlow
                dIst(1,13+numIn*5+numEx*6+i*5) = (x(1)-x(2)-thisNeuron.E_NMDA)*thisNeuron.g_ExSlowReceptor*B;
                dIst(2,13+numIn*5+numEx*6+i*5) = -dIst(1,13+numIn*5+numEx*6+i*5);
            end
        end
        
        
        
        function numMCs = determineNumberOfMCs(thisNeuron)
            numKChannels = thisNeuron.numKChannels;
            
            numNaChannels = thisNeuron.numNaChannels;
            
            numInReceptors = thisNeuron.numInReceptors;
            numExReceptors = thisNeuron.numExReceptors;
            numExSlowReceptors = thisNeuron.numExSlowReceptors;
            
            numMCs = numKChannels + numNaChannels + numInReceptors + numExReceptors + numExSlowReceptors;
        end
        
        
        
        function registerMCEffects(thisNeuronCh, aChannel, oldStateVector, newStateVector)
            if isa(aChannel,'KChannelMC') == true
                
                if any(oldStateVector ~= newStateVector)
                    thisNeuronCh.totalKConductance = newStateVector(5,1)*aChannel.g_K;
                end
                
            elseif isa(aChannel,'NaChannelMC') == true
                
                if any(oldStateVector ~= newStateVector)
                    thisNeuronCh.totalNaConductance = newStateVector(8,1)*aChannel.g_Na;
                end
                
            elseif isa(aChannel,'inReceptorMC') == true
                
                if any(oldStateVector ~= newStateVector)
                    thisNeuronCh.totalInReceptorConductance ...
                        = thisNeuronCh.totalInReceptorConductance + sum(newStateVector(4:5,1)-oldStateVector(4:5,1))*aChannel.g_InReceptor;
                end
                
            elseif isa(aChannel,'exReceptorMC') == true
                
                if any(oldStateVector ~= newStateVector)
                    thisNeuronCh.totalExReceptorConductance ...
                        = thisNeuronCh.totalExReceptorConductance + (newStateVector(6,1)-oldStateVector(6,1))*aChannel.g_ExReceptor;
                end
                
            elseif isa(aChannel,'exSlowReceptorMC') == true
                
                if any(oldStateVector ~= newStateVector)
                    thisNeuronCh.totalExSlowReceptorConductance ...
                        = thisNeuronCh.totalExSlowReceptorConductance + (newStateVector(5,1)-oldStateVector(5,1))*aChannel.g_ExSlowReceptor;
                end
                
            else
                error('This channel is not a valid channel object!');
            end
        end
        
        
        function synapse(thisNeuron, sinkNeuron, varargin)
            
            if isa(sinkNeuron, 'neuronMC') == true
                for j=1:1:sinkNeuron.numMCs
                    
                    if (isa(sinkNeuron.MCs{j}, varargin{:}) == true) && (isempty(sinkNeuron.MCs{j}.sourceNeuron) == true)
                        sinkNeuron.MCs{j}.sourceNeuron = thisNeuron;
                        break;
                    end
                end
                
            elseif isa(sinkNeuron, 'neuronHH') == true
                if strcmp(varargin{:}, 'exReceptor')
                    for i = 1:1:sinkNeuron.numExReceptors
                        if isempty(sinkNeuron.exSourceNeuron{1,i}) == true
                            sinkNeuron.exSourceNeuron{1,i} = thisNeuron;
                            break;
                        end
                    end
                    
                elseif strcmp(varargin{:}, 'exSlowReceptor')
                    for i = 1:1:sinkNeuron.numExSlowReceptors
                        if isempty(sinkNeuron.exSlowSourceNeuron{1,i}) == true
                            sinkNeuron.exSlowSourceNeuron{1,i} = thisNeuron;
                            break;
                        end
                    end
                    
                elseif strcmp(varargin{:}, 'inReceptor')
                    for i = 1:1:sinkNeuron.numInReceptors
                        if isempty(sinkNeuron.inSourceNeuron{1,i}) == true
                            sinkNeuron.inSourceNeuron{1,i} = thisNeuron;
                            break;
                        end
                    end
                end
            end
        end
        
        
        function [dRateVectorCross] = synapseRates(neuron1, neuron2, receptor)
            
            if isa(neuron2, 'neuronMC') == true
                sourceNeuron = neuron1;
                sinkNeuron = neuron2;
            elseif isa(neuron2, 'neuronHH') == true
                sourceNeuron = neuron2;
                sinkNeuron = neuron1;
            end
            
            numTransitionsSink = receptor.numTransitions;
            dRateVectorCross = zeros(numTransitionsSink, 2);
            
            x(1)=sourceNeuron.Xi(1);
            x(2)=sourceNeuron.Xi(2);
            
            Vpre=x(1)-x(2);
            
            if isa(receptor, 'inReceptorMC') == true
                
                Rb1=20*10^6;
                Rb2=10*10^6;
                
                dRateVectorCross(1,1) = Rb1 * ((-2.84*10^-3)/(1+exp(-(Vpre*10^3-2)/5))^2)*(-200)*exp(-(Vpre*10^3-2)/5);
                dRateVectorCross(3,1) = Rb2 * ((-2.84*10^-3)/(1+exp(-(Vpre*10^3-2)/5))^2)*(-200)*exp(-(Vpre*10^3-2)/5);
                dRateVectorCross(1,2) = -dRateVectorCross(1,1);
                dRateVectorCross(3,2) = -dRateVectorCross(3,1);
                
            elseif isa(receptor, 'exReceptorMC') == true
                
                Rb=13*10^6;
                
                dRateVectorCross(1,1) = Rb * ((-2.84*10^-3)/(1+exp(-(Vpre*10^3-2)/5))^2)*(-200)*exp(-(Vpre*10^3-2)/5);
                dRateVectorCross(3,1) = Rb * ((-2.84*10^-3)/(1+exp(-(Vpre*10^3-2)/5))^2)*(-200)*exp(-(Vpre*10^3-2)/5);
                dRateVectorCross(1,2) = -dRateVectorCross(1,1);
                dRateVectorCross(3,2) = -dRateVectorCross(3,1);

            elseif isa(receptor, 'exSlowReceptorMC') == true
                
                Rb=5*10^6*0.2;
                
                dRateVectorCross(1,1) = Rb * ((-2.84*10^-3)/(1+exp(-(Vpre*10^3-2)/5))^2)*(-200)*exp(-(Vpre*10^3-2)/5);
                dRateVectorCross(3,1) = Rb * ((-2.84*10^-3)/(1+exp(-(Vpre*10^3-2)/5))^2)*(-200)*exp(-(Vpre*10^3-2)/5);
                dRateVectorCross(1,2) = -dRateVectorCross(1,1);
                dRateVectorCross(3,2) = -dRateVectorCross(3,1);
            end
            
        end
        
    end
    
    methods (Access = protected)
        
        
        function [Q, I, J, dQ, dI] = eqEval(thisNeuronCh, x, t)
            
            C_m = thisNeuronCh.c_m*thisNeuronCh.Area; %memrane capacitance
            q = C_m*(x(1)-x(2));
            
            Q = [ q;
                -q;];
            
            dQ = [C_m, -C_m;
                -C_m, C_m];
            
            G_K  = thisNeuronCh.totalKConductance;  %total K channel conductance
            G_Na = thisNeuronCh.totalNaConductance; %total Na channel conductance
            G_In  = thisNeuronCh.totalInReceptorConductance; %total inhibitory synapse conductance
            G_Ex = thisNeuronCh.totalExReceptorConductance;  %total excitatory synapse conductance
            G_ExSlow = thisNeuronCh.totalExSlowReceptorConductance;  %total excitatory slow synapse conductance
            G_L  = thisNeuronCh.g_L*thisNeuronCh.Area;  %leakage resistance
            %B = 1;% 1/(1+exp(-0.062*1e3*(x(1)-x(2)))*(1.5/3.57));
            %dB = 0;%((1+exp(-0.062*1e3*(x(1)-x(2)))*(1.5/3.57))^-2)*(0.062*1e3*exp(-0.062*1e3*(x(1)-x(2)))*(1.5/3.57));
            
            I(1,1) = ((x(1)-x(2)-thisNeuronCh.E_K)*G_K) + ((x(1)-x(2)-thisNeuronCh.E_Na)*G_Na) + ((x(1)-x(2)-thisNeuronCh.E_L)*G_L) +...
                ((x(1)-x(2)-thisNeuronCh.E_Cl)*G_In) + ((x(1)-x(2)-thisNeuronCh.E_AMPA)*G_Ex) + ((x(1)-x(2)-thisNeuronCh.E_NMDA)*G_ExSlow);
            I(2,1) = -I(1,1);
            
            dI(1,1) = G_K + G_Na + G_L + G_In + G_Ex + G_ExSlow;
            dI(1,2) = -dI(1,1);
            
            dI(2,1) = -dI(1,1);
            dI(2,2) = -dI(1,2);
            
            J = zeros(2,1);
        end
        
        function newNeuronCh = copyElement(thisNeuron)
            % this method is necessary to avoid creating new neurons
            % with identical channel constellations when copying them. Ideally
            % this would go under circuitComponentMC but then we'd
            % have two differing definitions of the same method!
            Area = thisNeuron.Area;
            numExInputs = thisNeuron.numExReceptors;
            numExSlowInputs = thisNeuron.numExSlowReceptors;
            numInInputs = thisNeuron.numInReceptors;
            newNeuronCh = neuronMC(Area, numExInputs, numExSlowInputs, numInInputs, thisNeuron.name);
            
        end
        
        
        
        function initMCs(thisNeuron)
            %KStateVector = floor(100*[36; 36; 36.; 36; 36]);
            %KStateVector = 0.1*[58.615; 75.914; 36.869; 11.937; 00.644];
            KStateVector = [4.0176; 7.4382; 5.1641; 1.5934; 0.1843];
            %KStateVector = [12306.4894402291;47568.6627862060;68950.6240383112;44419.3290721907;10730.8946630640];
            KStateChangeMatrix = [-1 1 0 0 0 0 0 0;
                                   1 -1 -1 1 0 0 0 0;
                                   0 0 1 -1 -1 1 0 0;
                                   0 0 0 0 1 -1 -1 1;
                                   0 0 0 0 0 0 1 -1];
            
            %NaStateVector = floor(100*[75; 75; 75; 75; 75; 75; 75; 75]);
            %NaStateVector = 0.1*[135.126; 12.064; 00.359; 00.004; 414.344; 36.992; 01.101; 00.011];
            NaStateVector = [20.4580; 3.3947; 0.1877; 0.0034; 30.5934; 5.0765; 0.2807; 0.0051];
            %NaStateVector = [352458.823155750;54456.8425082606;2804.62856364244;48.1478774312694;163623.468942195;25280.7322797944;1302.00477285974;22.3519000695330];
            NaStateChangeMatrix = [-1 1 0 0 0 0 0 0 0 0 0 0 -1 1 0 0 0 0 0 0;
                                    1 -1 -1 1 0 0 0 0 0 0 0 0 0 0 -1 1 0 0 0 0;
                                    0 0 1 -1 -1 1 0 0 0 0 0 0 0 0 0 0 -1 1 0 0;
                                    0 0 0  0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 -1 1;
                                    0 0 0  0 0 0 -1 1 0 0 0 0 1 -1 0 0 0 0 0 0;
                                    0 0 0  0 0 0 1 -1 -1 1 0 0 0 0 1 -1 0 0 0 0;
                                    0 0 0  0 0 0 0 0 1 -1 -1 1 0 0 0 0 1 -1 0 0;
                                    0 0 0  0 0 0 0 0 0 0 1 -1 0 0 0 0 0 0 1 -1];
            
            %inReceptorStateVector = 0.1*[60; 80; 60; 00; 00];
            inReceptorStateVector = [19.9995; 0.0003; 0.00001; 0.0001; 0.00001]/20;
            inReceptorStateChangeMatrix = [-1 1 0 0 0 0 0 0;
                                            1 -1 -1 1 -1 1 0 0;
                                            0 0 1 -1 0 0 -1 1;
                                            0 0 0 0 1 -1 0 0;
                                            0 0 0 0 0 0 1 -1];
            
            %exReceptorStateVector = 0.1*[200; 0; 0; 0; 0; 0];
            exReceptorStateVector = [1.2631; 1.2164; 17.5201; 0.00001; 0.00001; 0.00001]/20;
            exReceptorStateChangeMatrix = [-1 1 0 0 0 0 0 0 0 0;
                                            1 -1 -1 1 0 0 -1 1 0 0;
                                            0 0 0 0 0 0 1 -1 0 0;
                                            0 0 1 -1 -1 1 0 0 -1 1;
                                            0 0 0 0 0 0 0 0 1 -1;
                                            0 0 0 0 1 -1 0 0 0 0];
            %exSlowReceptorStateVector = 0.15*[60; 60; 60; 20; 0];
            exSlowReceptorStateVector = [27.0831; 2.6590; 0.00005; 0.2576; 0.00015]/30;
            exSlowReceptorStateChangeMatrix = [-1 1 0 0 0 0 0 0;
                                                1 -1 -1 1 0 0 0 0;
                                                0 0 1 -1 -1 1 -1 1;
                                                0 0 0 0 1 -1 0 0;
                                                0 0 0 0 0 0 1 -1];
            
            N = thisNeuron.determineNumberOfMCs();
            N_K = thisNeuron.numKChannels;
            N_Na = thisNeuron.numNaChannels;
            N_In = thisNeuron.numInReceptors;
            N_Ex = thisNeuron.numExReceptors;
            N_ExSlow = thisNeuron.numExSlowReceptors;
            
            n=0;
            m=N_K;
            for i = n+1:1:m
                thisNeuron.numMCs = i;
                ch = KChannelMC(KStateVector, KStateChangeMatrix, thisNeuron);
                thisNeuron.MCs{i} = ch;
                ch.setNumber(i);
            end
            thisNeuron.numStates = thisNeuron.numStates + N_K*ch.numStates;
            thisNeuron.numTransitions = thisNeuron.numTransitions + N_K*ch.numTransitions;
            thisNeuron.totalKConductance = KStateVector(5,1)*thisNeuron.g_K;
            
            
            n=n+N_K;
            m=m+N_Na;
            for i = n+1:1:m
                thisNeuron.numMCs = i;
                ch = NaChannelMC(NaStateVector, NaStateChangeMatrix, thisNeuron);
                thisNeuron.MCs{i} = ch;
                ch.setNumber(i);
            end
            thisNeuron.numStates = thisNeuron.numStates + N_Na*ch.numStates;
            thisNeuron.numTransitions = thisNeuron.numTransitions + N_Na*ch.numTransitions;
            thisNeuron.totalNaConductance = NaStateVector(8,1)*thisNeuron.g_Na;
            
            
            n=n+N_Na;
            m=m+N_In;
            for i = n+1:1:m
                thisNeuron.numMCs = i;
                ch = inReceptorMC(inReceptorStateVector, inReceptorStateChangeMatrix, thisNeuron);
                thisNeuron.MCs{i} = ch;
                ch.setNumber(i);
            end
            thisNeuron.numStates = thisNeuron.numStates + N_In*ch.numStates;
            thisNeuron.numTransitions = thisNeuron.numTransitions + N_In*ch.numTransitions;
            thisNeuron.totalInReceptorConductance = N_In*(inReceptorStateVector(4,1) + inReceptorStateVector(5,1))*thisNeuron.g_InReceptor;
            
            
            n=n+N_In;
            m=m+N_Ex;
            for i = n+1:1:m
                thisNeuron.numMCs = i;
                ch = exReceptorMC(exReceptorStateVector, exReceptorStateChangeMatrix, thisNeuron);
                thisNeuron.MCs{i} = ch;
                ch.setNumber(i);
            end
            thisNeuron.numStates = thisNeuron.numStates + N_Ex*ch.numStates;
            thisNeuron.numTransitions = thisNeuron.numTransitions + N_Ex*ch.numTransitions;
            thisNeuron.totalExReceptorConductance = N_Ex*exReceptorStateVector(6,1)*thisNeuron.g_ExReceptor;
            
            
            n=n+N_Ex;
            m=m+N_ExSlow;
            for i = n+1:1:m
                thisNeuron.numMCs = i;
                ch = exSlowReceptorMC(exSlowReceptorStateVector, exSlowReceptorStateChangeMatrix, thisNeuron);
                thisNeuron.MCs{i} = ch;
                ch.setNumber(i);
            end
            thisNeuron.numStates = thisNeuron.numStates + N_ExSlow*ch.numStates;
            thisNeuron.numTransitions = thisNeuron.numTransitions + N_ExSlow*ch.numTransitions;
            thisNeuron.totalExSlowReceptorConductance = N_ExSlow*exSlowReceptorStateVector(5,1)*thisNeuron.g_ExSlowReceptor;
        end
        
        
    end
    
end
