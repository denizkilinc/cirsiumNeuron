classdef transientSolverMC < handle
    % TRANSIENTSOLVERMC time-domain Monte Carlo solver for a neuronal 
    % circuit. This solver is can used for both MC-based and SDE-based 
    % stochastic Monte Carlo type simulations of a neuronal circuit.
    
    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------
    
    properties (Constant)
        TEROH = 0.1; % MC event rate overhead
        
        keyNames = {'solverType'
            'linSolverType'
            'currentVariables'
            't0'
            'tend'
            'timeStep'
            'breakPoints'
            'reltol'
            'abstol_v'
            'abstol_c'
            'abstol_q'
            'maxOrder'
            'maxStep'
            'chargeScaleFactor'
            'y0'
            'yp0'
            's0'
            'sp0'
            'seq0'
            'showSolution'
            'lmax'};
        validSolverTypes = {'IDA-SSA', 'SEQ-TRPZ'};
        validLinSolverTypes = {'Dense', 'GMRES', 'BiCGStab'};
        validOnOff = {'on', 'off'};
    end
    
    
    properties (SetAccess = private)
        circuit = circuit.empty;
        solverType = 'IDA-SSA';
        linSolverType = 'GMRES';
        currentVariables = 'off';
        showSolution = 'on';
        resFunIC;
        jacFunIC;
        precSolveFunIC;
        resFun;
        jacFun;
        precSolveFun;
        QFun;
        IFun;
        JFun;
        dQFun;
        dIFun;
        dQstFun;
        dIstFun;
        stampFun;
        ttcFun;
        t0 = 0;
        breakPoints;
        reltol = 1e-3;
        abstol_v = 1e-6;
        abstol_c = 1e-9;
        abstol_q = 1e-21;
        restol = 1e-6;
        tolVecIC;
        maxOrder = 5;
        maxStep = inf;
        timeStep = 1e-12;
        chargeScaleFactor = 1e9;
        y0;
        yp0;
        y0I;
        yp0I;
        s0;
        sp0;
        seq0;
        seqp0;
        ynow;
        ypnow;
        Qprev;
        tprev;
        tnow;
        tend;
        tolVec;
        IDAOptions;
        Y;
        T;
        SC;
        SV;
        
        TE;
        TC;
        PROPENSITYVECTOR;
        RATEVECTOR;
        
        stats;
        numVars;
        numVoltageVars;
        numNewtonIters;
        condNumbers;
        nextStopTime;
        nextBreakPoint;
        
        numMCs;
        
        numStates;
        numTransitions;
        
        snow;
        lmax;
        lnow;
        
        propensityVectorNow;
        rateVectorNow;
        
        loverhead = transientSolverMC.TEROH;
        S;
        Seq;
        eventTimes;
        nextEventTime;
        
        firedTraps;
        firedMCs;
        
        numBreaks = 0;
        numIters = 0;
        initialDerivative;
        calcCurrents = false;
        lastState;
        
    end
    
    methods
        
        function obj = transientSolverMC(circ, varargin)
            if nargin > 0
                if nargin > 1
                    obj.setOptions(varargin{:});
                end
                obj.setCircuit(circ);
            end
        end
        
        function setOptions(thisSolver, varargin)
            if rem(nargin-1, 2) ~= 0
                error('Options structure is missing a key or a value!');
            end
            for i=1:2:nargin-1
                key = varargin{i};
                if ~ischar(key)
                    error('Options must consist of key-value pairs');
                end
                val = varargin{i + 1};
                try
                    thisSolver.(key) = val;
                catch e
                    disp(e);
                    error('%s is an unrecognized key!', key);
                end
            end
            
            % check if the options are valid
            
            if all(strcmp(thisSolver.solverType, transientSolverMC.validSolverTypes) == 0)
                msgStr = sprintf('Unrecognized solver type: %s', thisSolver.solverType);
                msgStr = sprintf('%s\nValid choices are:', msgStr);
                for i=1:numel(transientSolverMC.validSolverTypes)
                    msgStr = sprintf('%s\n%d. %s', msgStr, i, transientSolverMC.validSolverTypes{i});
                end
                error(msgStr);
            end
            
            if all(strcmp(thisSolver.linSolverType, transientSolverMC.validLinSolverTypes) == 0)
                msgStr = sprintf('Unrecognized linear solver type: %s', thisSolver.linSolverType);
                msgStr = sprintf('%s\nValid choices are:', msgStr);
                for i=1:numel(transientSolverMC.validLinSolverTypes)
                    msgStr = sprintf('%s\n%d. %s', msgStr, i, transientSolverMC.validLinSolverTypes{i});
                end
                error(msgStr);
            end
            
            if all(strcmp(thisSolver.currentVariables, transientSolverMC.validOnOff) == 0)
                msgStr = sprintf('Unrecognized current variables option: %s', thisSolver.currentVariables);
                msgStr = sprintf('%s\nValid choices are:', msgStr);
                for i=1:numel(transientSolverMC.validOnOff)
                    msgStr = sprintf('%s\n%d. %s', msgStr, i, transientSolverMC.validOnOff{i});
                end
                error(msgStr);
            end
            
            if all(strcmp(thisSolver.showSolution, transientSolverMC.validOnOff) == 0)
                msgStr = sprintf('Unrecognized plotting option: %s', thisSolver.showSolution);
                msgStr = sprintf('%s\nValid choices are:', msgStr);
                for i=1:numel(transientSolverMC.validOnOff)
                    msgStr = sprintf('%s\n%d. %s', msgStr, i, transientSolverMC.validOnOff{i});
                end
                error(msgStr);
            end
            
            % did we get an initial derivative vector?
            thisSolver.initialDerivative = ~isempty(thisSolver.yp0);
            
            % if so we can calculate a steady state solution etc.
            
            % set the solver functions for initial condition calculations
            % with IDA
            resFunName = 'IDAResFun';
            precSolveFunName = 'IDAPrecSolveFun';
            switch thisSolver.linSolverType
                case 'GMRES'
                    jacFunName = 'IDAJacTimesVecFun';
                case 'BiCGStab'
                    jacFunName = 'IDAJacTimesVecFun';
                case 'Dense'
                    jacFunName = 'IDADJacFun';
            end
            
            if strcmp(thisSolver.currentVariables, 'off')
                thisSolver.calcCurrents = true;
            end
            
            if strncmp(thisSolver.solverType, {'SEQ'}, 3)
                resFunName = [resFunName, 'EMS'];
                jacFunName = [jacFunName, 'EMS'];
                precSolveFunName = [precSolveFunName, 'EMS'];
            end
            
            thisSolver.resFunIC = eval(['@thisSolver.' resFunName]);
            thisSolver.jacFunIC = eval(['@thisSolver.' jacFunName]);
            thisSolver.precSolveFunIC = eval(['@thisSolver.' precSolveFunName]);
            
            % set the solver functions
            
            switch thisSolver.solverType
                case {'IDA-SSA'}
                    thisSolver.resFun = thisSolver.resFunIC;
                    thisSolver.jacFun = thisSolver.jacFunIC;
                    thisSolver.precSolveFun = thisSolver.precSolveFunIC;
                case 'SEQ-TRPZ'
                    thisSolver.resFun = @thisSolver.SEQ_TRPZResFun;
                    thisSolver.jacFun = @thisSolver.SEQ_TRPZJacFun;
            end
            
        end
        
        function setCircuit(thisSolver, circ)
            if isa(circ, 'circuit') || isa(circ, 'subcircuit')
                thisSolver.circuit = circ;
            else
                error('The solver only accepts objects of type circuit or subcircuit');
            end
            
            % are the source currents included as variables?
            switch thisSolver.currentVariables
                case 'on'
                    thisSolver.numVars = circ.numVars;
                    thisSolver.numVoltageVars = circ.numIndepVoltVars;
                case 'off'
                    thisSolver.numVars = circ.numVarsnc;
                    thisSolver.numVoltageVars = thisSolver.numVars;
            end
            
            % set the number of MCs
            thisSolver.numMCs = circ.numMCs;
            thisSolver.numStates = circ.numStates;
            thisSolver.numTransitions = circ.numTransitions;
            
            if thisSolver.initialDerivative && isempty(thisSolver.sp0)
                thisSolver.sp0 = zeros(thisSolver.numStates,1);
            end
            
            % if there are MCs require a lambda_max
            if thisSolver.numMCs >0 && ...
                    strcmp(thisSolver.solverType, 'IDA-SSA') == true &&...
                    (isempty(thisSolver.lmax) || thisSolver.lmax == 0)
                error('Please specify a valid lmax!');
            end
            
            % check if initial conditions have the correct size
            
            if ~isempty(thisSolver.s0) && (thisSolver.numStates ~= length(thisSolver.s0))
                error('Size of s0 is not consistent with the number of traps!');
            end
            
            if ~isempty(thisSolver.seq0) && (thisSolver.numStates ~= length(thisSolver.seq0))
                error('Size of seq0 is not consistent with the number of traps!');
            end
            
            if isempty(thisSolver.y0)
                error(['please supply an initial condition vector or  an'...
                    'initial guess for the steady state solution']);
            elseif length(thisSolver.y0) ~= thisSolver.numVars
                error('Size of y0 is not consistent with the problem size!');
            end
            
            if (thisSolver.initialDerivative) && (length(thisSolver.yp0) ~= thisSolver.numVars)
                error('Size of yp0 is not consistent with the problem size!');
            end
            
            % construct tolerance vectors
            nV = thisSolver.numVars;
            nVV = thisSolver.numVoltageVars;
            nMCs = thisSolver.numMCs;
            nStates = thisSolver.numStates;
            Kchrg = thisSolver.chargeScaleFactor;
            
            thisSolver.tolVecIC = [ones(nVV,1).*thisSolver.abstol_v;
                ones(nV-nVV,1).*thisSolver.abstol_c;
                ones(nV,1)*thisSolver.abstol_q*thisSolver.chargeScaleFactor];
            
            switch thisSolver.solverType(1:3)
                case 'IDA'
                    thisSolver.tolVecIC = [ones(nVV,1).*thisSolver.abstol_v;
                        ones(nV-nVV,1).*thisSolver.abstol_c;
                        ones(nV,1)*thisSolver.abstol_q*Kchrg];
                    thisSolver.tolVec = thisSolver.tolVecIC;
                case 'SEQ'
                    thisSolver.tolVecIC = [ones(nVV,1).*thisSolver.abstol_v;
                        ones(nV-nVV,1).*thisSolver.abstol_c;
                        ones(nV,1)*thisSolver.abstol_q*Kchrg;
                        ones(nStates,1)*thisSolver.abstol_v];
                    thisSolver.tolVec = [ones(nVV,1).*thisSolver.abstol_v;
                        ones(nV-nVV,1).*thisSolver.abstol_c;
                        ones(2*nStates,1)*thisSolver.abstol_v];
                    thisSolver.restol = [ones(nVV,1)*thisSolver.abstol_q*Kchrg;
                        ones(nV-nVV+nStates,1)*thisSolver.abstol_v;
                        ones(nStates,1)*thisSolver.abstol_v];
            end
            
            % set native residual and jacobian functions
            switch thisSolver.currentVariables
                case 'on'
                    thisSolver.QFun = @circ.Q;
                    thisSolver.IFun = @circ.I;
                    thisSolver.JFun = @circ.J;
                    thisSolver.dQFun = @circ.dQ;
                    thisSolver.dIFun = @circ.dI;
                    thisSolver.dQstFun = @circ.dQst;
                    thisSolver.dIstFun = @circ.dIst;
                    thisSolver.stampFun = @circ.stamp;
                    if nMCs > 0
                        thisSolver.ttcFun = @circ.ttc;
                    end
                case 'off'
                    thisSolver.QFun = @circ.Qnc;
                    thisSolver.IFun = @circ.Inc;
                    thisSolver.JFun = @circ.Jnc;
                    thisSolver.dQFun = @circ.dQnc;
                    thisSolver.dIFun = @circ.dInc;
                    thisSolver.dQstFun = @circ.dQstnc;
                    thisSolver.dIstFun = @circ.dIstnc;
                    thisSolver.stampFun = @circ.stampnc;
                    if nMCs > 0
                        thisSolver.ttcFun = @circ.ttcnc;
                    end
            end
            
            
        end
        
        function setIDAOptions(thisSolver)
            nV = thisSolver.numVars;
            nStates = thisSolver.numStates;
            
            diffInd = [zeros(nV,1); ones(nV,1)];
            
            if strncmp(thisSolver.solverType, {'SEQ'}, 3)
                diffInd = [diffInd; ones(nStates,1)];
            end
            
            % set the end of times as the last break point
            thisSolver.breakPoints = [thisSolver.breakPoints,...
                thisSolver.tend];
            % set the first break point as the next break point
            thisSolver.nextBreakPoint = thisSolver.breakPoints(1);
            
            % set the solver functions in the IDAOptions structure
            thisSolver.IDAOptions = IDASetOptions(...
                'LinearSolver', thisSolver.linSolverType,...
                'PrecModule', 'UserDefined',...
                'VariableTypes', diffInd,...
                'RelTol', thisSolver.reltol,...
                'AbsTol', thisSolver.tolVecIC,...
                'MaxNumSteps', 500,...
                'MaxStep', thisSolver.maxStep,...
                'MaxOrder', thisSolver.maxOrder,...
                'JacobianFn', thisSolver.jacFunIC,...
                'PrecSolveFn', thisSolver.precSolveFunIC,...
                'StopTime', thisSolver.nextBreakPoint);
            
        end
        
        
        function solve(thisSolver)
            % this function calls the main solver loop.
            if isempty(thisSolver.circuit)
                error(['This solver does not have a circuit to solve! ',...
                    'Please add a circuit first.']);
            elseif thisSolver.circuit.isSealed == false
                error('The circuit is not sealed! Please seal it first.');
            end
            
            thisSolver.setIDAOptions();
            thisSolver.calcIC();
            
            switch thisSolver.solverType(1:3)
                case 'IDA'
                    thisSolver.solveIDA;
                case 'SEQ'
                    thisSolver.solveSEQ;
            end
            
            if strcmp(thisSolver.showSolution, 'on')
                thisSolver.displaySolution;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %           Consistent initial conditions calculation             %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function calcIC(thisSolver)
            y0 = thisSolver.y0;
            yp0 = thisSolver.yp0;
            t0 = thisSolver.t0;
            sp0 = thisSolver.sp0;
            Kchrg = thisSolver.chargeScaleFactor;
            
            % first set the initial trap states
            if thisSolver.numMCs > 0
                if isempty(thisSolver.s0)
                    thisSolver.s0 = zeros(thisSolver.numStates,1);
                    
                    k=0;
                    for i=1:1:thisSolver.numMCs
                        m = thisSolver.circuit.MCs{i}.numStates;
                        thisSolver.s0(k+1:k+m,1) = thisSolver.circuit.MCs{i}.stateVector;
                        k=k+m;
                    end
                end
                thisSolver.circuit.setMCStates(thisSolver.s0);
                s0 = thisSolver.s0;
            end
            
            y02 = [y0; Kchrg*thisSolver.QFun(y0,t0)];
            if thisSolver.initialDerivative
                % user wants the solution with specified derivatives,
                % e.g. steady state solution
                yp02 = [yp0; Kchrg*(thisSolver.dQFun(y0,t0)*yp0)];
            else
                % we are only given the algebraic initial conditions
                yp02 = [zeros(thisSolver.numVars,1);
                    Kchrg*(thisSolver.JFun(y0,t0) - thisSolver.IFun(y0,t0))];
            end
            
            if strncmp(thisSolver.solverType, {'SEQ'}, 3)
                y02 = [y02; s0];
                if thisSolver.initialDerivative
                    % user wants the solution with spcified derivatives,
                    % e.g. steady state solution
                    yp02 = [yp02; sp0];
                else
                    % we are only given the algebraic initial conditions
                    yp02 = [yp02;
                        zeros(thisSolver.numStates,1)];
                end
            end
            y02 = full(y02);
            yp02 = full(yp02);
            
            IDAInit(thisSolver.resFunIC, t0, y02, yp02, thisSolver.IDAOptions);
            
            t1 = thisSolver.t0 + thisSolver.timeStep;
            
            tic
            if thisSolver.initialDerivative
                [~, thisSolver.y0I, thisSolver.yp0I] = IDACalcIC(t1, 'FindAll');
            else
                [~, thisSolver.y0I, thisSolver.yp0I] = IDACalcIC(t1, 'FindAlgebraic');
            end
            
            
            fprintf(1,'%g seconds for consistent initial conditions.\n', toc);
            
            % Free the IDA memory if we will not use it later
            if ~strncmp(thisSolver.solverType, 'IDA', 3)
                IDAFree;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %                    IDA - SSA solution                           %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function solveIDA(thisSolver)
            if strcmp(thisSolver.solverType, 'IDA-SSA')
                ssa = true;
            else
                ssa = false;
            end
            
            circ = thisSolver.circuit;
            
            % equation sizes
            nV = thisSolver.numVars;
            nTr = thisSolver.numMCs;
            nStates = thisSolver.numStates;
            
            bufferSize = 5e5;
            matrixSizeCheckTime = bufferSize;
            matrixSizeCheckFlip = bufferSize;
            
            % initial conditions
            y0I = thisSolver.y0I;
            yp0I = thisSolver.yp0I;
            t0 = thisSolver.t0;
            s0 = thisSolver.s0;
            tend = thisSolver.tend;
            
            % current state of the solver
            thisSolver.tnow = t0;
            thisSolver.ynow = y0I;
            thisSolver.ypnow = yp0I;
            thisSolver.snow = s0;
            
            % initialization of event times
            thisSolver.nextEventTime = Inf;
            thisSolver.numBreaks = 0;
            
            % initialization of solver outputs
            if ssa == true && nTr == 0
                thisSolver.Y = zeros(2*nV, bufferSize);
                thisSolver.T = zeros(1, bufferSize);
            elseif ssa == true && nTr > 0
                thisSolver.Y = zeros(2*nV, bufferSize);
                thisSolver.T = zeros(1, bufferSize);
                thisSolver.S = zeros(nStates, bufferSize);
                thisSolver.S(:,1) = s0;
            elseif ssa == false
                thisSolver.Y = zeros(2*nV + nStates, bufferSize);
                thisSolver.T = zeros(1, bufferSize);
            end
            
            thisSolver.Y(:,1) = y0I;
            thisSolver.T(:,1) = t0;
            
            goBack = false;
            
            % initial event rates
            if nTr > 0
                [rateVectorNow] = thisSolver.ttcFun(y0I(1:nV), t0, s0);
                thisSolver.RATEVECTOR = zeros(circ.numTransitions, bufferSize);
                thisSolver.RATEVECTOR(:,1) = rateVectorNow;
            end
            
            %field_snow = 'snow';
            
            % initial resetting fall-backs
            thisSolver.lastState = struct('tnow', t0,...
                'ynow', y0I,...
                'ypnow', yp0I,...
                'snow', s0,...
                'numBreaks', 0,...
                'numIters', 0,...
                'rng', rng());
            
            % do we calculate currents externally?
            calcCurrents = thisSolver.calcCurrents;
            
            if calcCurrents
                thisSolver.SV = zeros(circ.numVoltageSources, bufferSize);
                thisSolver.SC = zeros(circ.numVoltageSources, bufferSize);
                
                [SVt, SCt] = circ.getSourceVariables(...
                    y0I(1:nV), t0, [], [],...
                    yp0I(nV+1:2*nV));
                thisSolver.SV(:,1) = SVt;
                thisSolver.SC(:,1) = SCt;
                thisSolver.tprev = t0;
                thisSolver.Qprev = thisSolver.circuit.getFullQ(y0I(1:nV), t0);
                thisSolver.lastState.Qprev = thisSolver.Qprev;
                thisSolver.lastState.tprev = t0;
            end
            
            % set first event time
            if nTr > 0
                % update rates is necessary because we have to check if the
                % initial lmax is big enough.
                goBack = thisSolver.updateEventRates();
                if ssa == true
                    thisSolver.eventTimes = zeros(1,bufferSize);
                    thisSolver.firedMCs = zeros(2,bufferSize);
                    
                    thisSolver.setNextEventTime();
                end
            end
            
            % loop internals
            thisSolver.numIters = 0;
            titer = thisSolver.t0;
            
            % begin calculating the solution
            if goBack == false
                thisSolver.displayTime();
            end
            
            timeInd=2;
            flipInd=1;
            while thisSolver.tnow < thisSolver.tend
                % go on advancing with single steps until we hit a
                % break point
                % advance one time step
                [~, t, y] = IDASolve(tend, 'OneStep');
                
                % get the derivative vector for this time point
                thisSolver.ypnow = IDAGet('DerivSolution', t, 1);
                
                % update solver outputs
                thisSolver.ynow = y;
                thisSolver.tnow = t;
                thisSolver.Y(:,timeInd) = y;
                thisSolver.T(:,timeInd) = t;
                
                if ssa == true && nTr > 0
                    thisSolver.RATEVECTOR(:,timeInd) = thisSolver.rateVectorNow;
                end
                
                % calculate currents externally if necessary and save them.
                if calcCurrents == true
                    [SVt, SCt] = circ.getSourceVariables(y(1:nV), t, ...
                        thisSolver.Qprev, ...
                        thisSolver.tprev);
                    thisSolver.SC(:,timeInd) = SCt;
                    thisSolver.SV(:,timeInd) = SVt;
                    thisSolver.Qprev = thisSolver.circuit.getFullQ(y(1:nV), t);
                    thisSolver.tprev = t;
                end
                
                % update loop internals and show time point to the user.
                thisSolver.numIters = thisSolver.numIters + 1;
                if (t - titer)/tend * 100 > 1 || thisSolver.numIters > 1000
                    thisSolver.displayTime('override');
                    titer = t;
                    thisSolver.numIters = 0;
                end
                
                % trap related updates
                if nTr > 0
                    % TODO: this resets to t=0 even if ssa==false
                    goBack = thisSolver.updateEventRates();
                    
                    if ssa == true
                        if goBack == false;
                            if thisSolver.tnow == thisSolver.nextEventTime
                                thisSolver.fireMCEvent(flipInd);
                                flipInd=flipInd+1;
                                if flipInd > matrixSizeCheckFlip
                                    thisSolver.eventTimes = [thisSolver.eventTimes zeros(1,bufferSize)];
                                    thisSolver.firedMCs = [thisSolver.firedMCs zeros(2,bufferSize)];
                                    matrixSizeCheckFlip = matrixSizeCheckFlip+bufferSize;
                                end
                            end
                            thisSolver.S(:,timeInd) = thisSolver.snow;
                        end
                    end
                end
                
                if thisSolver.tnow == thisSolver.nextBreakPoint;
                    thisSolver.IDAStepBreakPoint();
                end
                
                
                timeInd=timeInd+1;
                if timeInd > matrixSizeCheckTime
                    if (ssa == true && nTr == 0) || (ssa == false)
                        thisSolver.Y = [thisSolver.Y zeros(2*nV + nStates, bufferSize)];
                        thisSolver.T = [thisSolver.T zeros(1, bufferSize)];
                    elseif ssa == true && nTr > 0
                        thisSolver.Y = [thisSolver.Y zeros(2*nV, bufferSize)];
                        thisSolver.T = [thisSolver.T zeros(1, bufferSize)];
                        thisSolver.S = [thisSolver.S zeros(nStates, bufferSize)];
                        thisSolver.RATEVECTOR = [thisSolver.RATEVECTOR zeros(circ.numTransitions, bufferSize)];
                    end
                    
                    if calcCurrents == true
                        thisSolver.SV = [thisSolver.SV zeros(circ.numVoltageSources, bufferSize)];
                        thisSolver.SC = [thisSolver.SC zeros(circ.numVoltageSources, bufferSize)];
                    end
                    matrixSizeCheckTime = matrixSizeCheckTime+bufferSize;
                end
            end
            
            %exclude extra zeros
            if (ssa == true && nTr == 0) || (ssa == false)
                thisSolver.Y = thisSolver.Y(:,1:timeInd-1);
                thisSolver.T = thisSolver.T(:,1:timeInd-1);
            elseif ssa == true && nTr > 0
                thisSolver.Y = thisSolver.Y(:,1:timeInd-1);
                thisSolver.T = thisSolver.T(:,1:timeInd-1);
                thisSolver.S = thisSolver.S(:,1:timeInd-1);
                thisSolver.RATEVECTOR = thisSolver.RATEVECTOR(:,1:timeInd-1);
                thisSolver.firedMCs = thisSolver.firedMCs(:,1:flipInd-1);
                thisSolver.eventTimes = thisSolver.eventTimes(:,1:flipInd-1);
            end
            
            if calcCurrents == true
                thisSolver.SV = thisSolver.SV(:,1:timeInd-1);
                thisSolver.SC = thisSolver.SC(:,1:timeInd-1);
            end
            
            if ssa == false
                thisSolver.S = thisSolver.Y(2*nV+1:2*nV+nStates,:);
                thisSolver.Y = thisSolver.Y(1:2*nV,:);
            end
            
            % end solution save IDA state
            thisSolver.stats = IDAGetStats;
            IDAFree;
        end
        
        function goBack = updateEventRates(thisSolver)
            goBack = false;
            snow = thisSolver.snow;
            tnow = thisSolver.tnow;
            
            [rateVectorNow, dRateVectorNow, dRateVectorStateNow] = thisSolver.ttcFun(thisSolver.ynow(1:thisSolver.numVars), tnow, snow);
            
            [propensityVectorNow] = thisSolver.propensityFun(rateVectorNow, dRateVectorNow, dRateVectorStateNow, snow);
            
            lnow=propensityVectorNow;
            
            thisSolver.rateVectorNow = rateVectorNow;
            
            thisSolver.lnow = lnow; %here note that lnow keeps the total event rates of each MC
            
            if (sum(lnow) > thisSolver.lmax)
                thisSolver.lmax = sum(lnow)*(1+thisSolver.TEROH);
                thisSolver.resetToLastEventTime();
                goBack = true;
            end
        end
        
        function setNextEventTime(thisSolver)
            dt = thisSolver.randEventTime();
            thisSolver.nextEventTime = thisSolver.tnow + dt;
            thisSolver.IDASetNextStopTime();
        end
        
        function ret = randEventTime(thisSolver)
            ret = exprnd(1/(thisSolver.lmax * (1 + thisSolver.loverhead)));
        end
        
        function IDASetNextStopTime(thisSolver)
            if thisSolver.nextEventTime < thisSolver.nextBreakPoint
                thisSolver.nextStopTime = thisSolver.nextEventTime;
            else
                thisSolver.nextStopTime = thisSolver.nextBreakPoint;
            end
            thisSolver.IDAOptions.StopTime = thisSolver.nextStopTime;
            IDASet('StopTime', thisSolver.nextStopTime);
        end
        
        function fireMCEvent(thisSolver, flipInd)
            rn = rand()*thisSolver.lmax;
            ind = find(cumsum(thisSolver.lnow) > rn, 1, 'first');
            
            dummyEvent = isempty(ind);
            if dummyEvent == false
                dummyCheck = true;
                while(dummyCheck == true)
                    eventInd=0;
                    for i=1:1:thisSolver.numMCs
                        m = thisSolver.circuit.MCs{i}.numTransitions;
                        for j=1:1:m
                            eventInd = eventInd + 1;
                            if eventInd == ind
                                dummyCheck = false;
                                break;
                            end
                        end
                        if eventInd == ind
                            dummyCheck = false;
                            break;
                        end
                    end
                    break;
                end
                
                indMC=i;
                indEvent=j;
                
                thisSolver.flipMCState(indMC,indEvent);
                thisSolver.firedMCs(:,flipInd) = [indMC; indEvent];
                %fprintf(1,'\nIn MC %d transition %d occured at t = %g.\n', indMC, indEvent, thisSolver.tnow);
                lastState = thisSolver.lastState;
                lastState.tnow = thisSolver.tnow;
                lastState.ynow = thisSolver.ynow;
                lastState.ypnow = thisSolver.ypnow;
                lastState.numIters = thisSolver.numIters;
                lastState.numBreaks = thisSolver.numBreaks;
                
                if thisSolver.calcCurrents == true
                    lastState.Qprev = thisSolver.Qprev;
                    lastState.tprev = thisSolver.tprev;
                end
                
                thisSolver.lastState = lastState;
            else
                thisSolver.firedMCs(:,flipInd) = [0; 0];
                %fprintf(1,'\nDummy event at t = %g.\n', thisSolver.tnow);
            end
            
            thisSolver.eventTimes(:,flipInd) = thisSolver.tnow;
            thisSolver.setNextEventTime();
            
            %thisSolver.displayTime();
            % calculate the starting point after state change
            if dummyEvent == false
                IDAReInit(thisSolver.tnow, thisSolver.ynow, thisSolver.ypnow, ...
                    thisSolver.IDAOptions);
                [~, thisSolver.ynow, thisSolver.ypnow] = ...
                    IDACalcIC(thisSolver.tnow + thisSolver.timeStep, 'FindAlgebraic');
            end
        end
        
        
        function flipMCState(thisSolver, indMC, indEvent)
            
            thisSolver.circuit.flipMCState(indMC,indEvent);
            
            rateVectorNow = thisSolver.rateVectorNow;
            
            snow=[];
            for i=1:1:thisSolver.circuit.numMCs
                snow = [snow; thisSolver.circuit.MCs{i}.stateVector];
            end
            
            thisSolver.snow = snow;
        end
        
        
        function IDAStepBreakPoint(thisSolver)
            % we either hit a break point or the end of times
            if thisSolver.tnow < thisSolver.tend
                thisSolver.numBreaks = thisSolver.numBreaks + 1;
                
                fprintf(1,'\nHit a break point at t = %g.\n', thisSolver.nextBreakPoint);
                thisSolver.nextBreakPoint = thisSolver.breakPoints(thisSolver.numBreaks+1);
                thisSolver.IDASetNextStopTime();
                
                disp('Reinitializing...');
                
                IDAReInit(thisSolver.tnow, thisSolver.ynow, thisSolver.ypnow, thisSolver.IDAOptions);
                [~, thisSolver.ynow, thisSolver.ypnow] = ...
                    IDACalcIC(thisSolver.tnow+thisSolver.timeStep, 'FindAlgebraic');
                thisSolver.displayTime();
            end
        end
        
        
        function resetToLastEventTime(thisSolver)
            fprintf(1,'\nReset to last trap event time at t = %g.\n', thisSolver.tnow);
            fprintf(1,'new lmax: %g\n', thisSolver.lmax);
            
            lastState = thisSolver.lastState;
            ind = find(thisSolver.T == lastState.tnow);
            thisSolver.T = thisSolver.T(1:ind);
            thisSolver.Y = thisSolver.Y(:,1:ind);
            thisSolver.S = thisSolver.S(:,1:ind);
            
            thisSolver.RATEVECTOR = thisSolver.RATEVECTOR(:,1:ind);
            
            thisSolver.tnow = lastState.tnow;
            thisSolver.ynow = lastState.ynow;
            thisSolver.ypnow = lastState.ypnow;
            thisSolver.numIters = lastState.numIters;
            thisSolver.numBreaks = lastState.numBreaks;
            thisSolver.snow = lastState.snow;
            
            if thisSolver.calcCurrents == true
                thisSolver.Qprev = lastState.Qprev;
                thisSolver.tprev = lastState.tprev;
                thisSolver.SC = thisSolver.SC(:,1:ind);
                thisSolver.SV = thisSolver.SV(:,1:ind);
            end
            
            % reset the random number generator.
            rng(lastState.rng);
            
            thisSolver.nextBreakPoint = thisSolver.breakPoints(...
                find((thisSolver.breakPoints > lastState.tnow), 1));
            
            IDAReInit(thisSolver.tnow, thisSolver.ynow, thisSolver.ypnow, thisSolver.IDAOptions);
            [~, thisSolver.ynow, thisSolver.ypnow] = ...
                IDACalcIC(thisSolver.tnow+thisSolver.timeStep/1e2, 'FindAlgebraic');
            [thisSolver.rateVectorNow] = thisSolver.ttcFun(...
                thisSolver.ynow(1:thisSolver.numVars),...
                thisSolver.tnow, thisSolver.snow);
            thisSolver.setNextEventTime();
            thisSolver.displayTime();
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %           SDE Solution with Smooth Diffusion Algorithm          %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        function solveSEQ(thisSolver)
            % Langevin solution (predictor/corrector scheme
            % for trap states)
            circ = thisSolver.circuit;
            
            h = thisSolver.timeStep;
            nV = thisSolver.numVars;
            nTr = thisSolver.numMCs;
            nStates = thisSolver.numStates;
            nTransitions = thisSolver.numTransitions;
            stateChangeMatrix = circ.stateChangeMatrixi;
            stateChangeMatrixDiff = circ.stateChangeMatrixiDiff;
            
            y0I = thisSolver.y0I;
            yp0I = thisSolver.yp0I;
            seq0 = thisSolver.seq0;
            Kchrg = thisSolver.chargeScaleFactor;
            
            thisSolver.T = thisSolver.t0:h:thisSolver.tend;
            numT = numel(thisSolver.T);
            t = thisSolver.t0;
            thisSolver.tnow = thisSolver.t0;
            
            thisSolver.Y = zeros(nV,numT);
            v = y0I(1:nV);
            thisSolver.ynow = v;
            thisSolver.Y(:,1) = v;
            
            thisSolver.S = zeros(nStates, numT);
            s = y0I(2*nV+1:2*nV+nStates);
            thisSolver.snow = s;
            thisSolver.S(:,1) = s;
            
            thisSolver.Seq = zeros(nStates, numT);
            seq = seq0;
            thisSolver.Seq(:,1) = seq;
            
            thisSolver.PROPENSITYVECTOR = zeros(nTransitions, numT);
            
            [Qprev, Iprev, Jprev] = thisSolver.stampFun(v,t,s);            
            
            [rateVector, dRateVector, dRateVectorState] = ...
                thisSolver.ttcFun(v,t,s);
            
            [propensityVector] = ...
                thisSolver.propensityFun(rateVector, dRateVector, dRateVectorState, s);
            [propensityVectorEq] = ...
                thisSolver.propensityFun(rateVector, dRateVector, dRateVectorState, seq);
            
            thisSolver.numNewtonIters = zeros(1,numT);
            
            calcCurrents = thisSolver.calcCurrents;
            
            if calcCurrents
                thisSolver.SC = zeros(circ.numVoltageSources,numT);
                thisSolver.SV = zeros(circ.numVoltageSources,numT);
                [SVt, SCt] = circ.getSourceVariables(y0I(1:nV), t, [], [],...
                    yp0I(nV+1:2*nV));
                thisSolver.SC(:,1) = SCt;
                thisSolver.SV(:,1) = SVt;
            end
            
            thisSolver.displayTime();
            
            %first prediction for all variables
            dW = sqrt(h).*randn(nTransitions,1);
            %DW = sqrt(Lc + seq .* (Le - Lc)).*dW;
            
            DW = stateChangeMatrixDiff*(sqrt(propensityVectorEq).*dW);
            
            % predictor calculations
            sprdct = s + yp0I(2*nV+1:2*nV+nStates)*h + DW;
            vprdct = v + yp0I(1:nV)*h;
            %seqprdct = seq + (Lc - (Lc + Le).*seq)*h;
            seqprdct = seq + (stateChangeMatrix*propensityVectorEq)*h;
            
            % end predictor
            xprdct = [vprdct; sprdct; seqprdct];
            
            thisSolver.PROPENSITYVECTOR(:,1) = propensityVector;
            
            for i=2:numT
                if calcCurrents == true
                    QprevFull = thisSolver.circuit.getFullQ(v(1:nV), thisSolver.tnow);
                end
                
                t = thisSolver.T(i);
                thisSolver.tnow = t;
                
                try
                    [x, numIter] = thisSolver.newtSol(thisSolver.resFun(t, Qprev, Iprev, Jprev, h, [s;seq], propensityVector, propensityVectorEq, dW),...
                        thisSolver.jacFun(t, h), xprdct);
                catch exception
                    disp('Newton solver failed');
                    error(getReport(exception));
                    break;
                end
                
                v = x(1:nV);
                s = x(nV+1:nV+nStates);
                seq = x(nV+nStates+1:nV+2*nStates);
                
                % extract the charges and time constants before changing
                % MC states
                [Qprev, Iprev, Jprev] = thisSolver.stampFun(v, thisSolver.tnow);
                %[propensityVector] = thisSolver.ttcFun(v,t,s);
                %[propensityVectorEq] = thisSolver.ttcFun(v,t,seq);
                
                [rateVector, dRateVector, dRateVectorState] = ...
                    thisSolver.ttcFun(v,t,s);
                [propensityVector] = ...
                    thisSolver.propensityFun(rateVector, dRateVector, dRateVectorState, s);
                [propensityVectorEq] = ...
                    thisSolver.propensityFun(rateVector, dRateVector, dRateVectorState, seq);
                
                % save solver variables
                thisSolver.ynow = v;
                thisSolver.Y(:,i) = v;
                thisSolver.S(:,i) = s;
                thisSolver.Seq(:,i) = seq;
                thisSolver.PROPENSITYVECTOR(:,i) = propensityVector;
                thisSolver.numNewtonIters(i) = numIter;
                if calcCurrents == true
                    [SVt, SCt] = circ.getSourceVariables(v, thisSolver.tnow,...
                        QprevFull, thisSolver.T(i-1));
                    thisSolver.SC(:,i) = SCt;
                    thisSolver.SV(:,i) = SVt;
                end
                
                % sample the random increment in the SDE
                dW = sqrt(h).*randn(nTransitions,1);
                %DW = sqrt(Lc + seq .* (Le - Lc)).*dW;
                
                % end predictor
                xprdct = [v; s; seq];
                
                if floor(i/numT*100) > floor((i-1)/numT*100)
                    thisSolver.displayTime('override');
                end
            end
            
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %                    User Display functions                       %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function displaySolution(thisSolver)
            %plot results
            nV = thisSolver.numVars;
            nVV = thisSolver.numVoltageVars;
            figW = 562;
            figH = 504;
            scrnSize = get(0,'ScreenSize');
            scrnH = scrnSize(4);
            figPos = [150, scrnH-figH-50, figW, figH];
            
            figure('OuterPosition', figPos);
            plot(thisSolver.T,thisSolver.Y(1:nVV,:)); grid on;
            title('Membrane voltages');
            
            %switch thisSolver.currentVariables
            %    case 'on'
            %        figure('OuterPosition', figPos + [figW, 0, 0, 0]);
            %        plot(thisSolver.T,thisSolver.Y(nVV+1:nV,:)); grid on;
            %        title('Currents');
            %    case 'off'
            %        figure('OuterPosition', figPos + [figW, 0, 0, 0]);
            %        plot(thisSolver.T,thisSolver.SC); grid on;
            %        title('Source Currents');
            %end
            
            if thisSolver.numMCs > 0
                figure('OuterPosition', figPos + [figW, 0, 0, 0]);
                plot(thisSolver.T,thisSolver.S); grid on;
                title('Ion channel states');
                %if ~strncmp(thisSolver.solverType, 'IDA', 3)
                %    figure('OuterPosition', figPos + [figW, -figH, 0, 0]);
                %    plot(thisSolver.T,thisSolver.numNewtonIters); grid on;
                %    title('Newton iterations');
                %end
            end
        end
        
        function displayTime(thisSolver,varargin)
            if nargin > 1 && strcmp(varargin{1}, 'override')
                fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b\b\b%12e\n', thisSolver.tnow);
            else
                fprintf(1,'time: %12e\n', thisSolver.tnow);
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %   Residual and Jacobian functions for various solvers           %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [R, flag] = IDAResFun(thisSolver, t, y, yp)
            nV = thisSolver.numVars;
            
            Kchrg = thisSolver.chargeScaleFactor;
            
            x1 = y(1:nV);
            x2 = y(nV+1:end);
            %x1p = yp(1:N);
            x2p = yp(nV+1:end);
            
            [Q, I, J] = thisSolver.stampFun(x1,t);
            
            R = [(x2p/Kchrg + I - J);
                x2 - Q*Kchrg];
            flag = 0;
        end
        
        function [J, flag] = IDADJacFun(thisSolver, t, y, ~, ~, cj)
            nV = thisSolver.numVars;
            
            Kchrg = thisSolver.chargeScaleFactor;
            
            x1 = y(1:nV);
            %x2 = y(N+1:end);
            
            [~, ~, ~, dQ, dI] = thisSolver.stampFun(x1,t);
            
            J = full([ dI,      sparse(nV,nV);
                -dQ*Kchrg   speye(nV)] + cj*[sparse(nV,nV) speye(nV)/Kchrg;
                sparse(nV,nV) sparse(nV,nV)]);
            flag = 0;
        end
        
        function [Jv, flag] = IDAJacTimesVecFun(thisSolver, t, y, ~, ~, v, cj)
            nV = thisSolver.numVars;
            
            Kchrg = thisSolver.chargeScaleFactor;
            
            x1 = y(1:nV);
            %x2 = y(N+1:end);
            
            [~, ~, ~, dQ, dI] = thisSolver.stampFun(x1,t);
            
            J = [ dI sparse(nV,nV);
                -dQ*Kchrg   speye(nV)] + cj*[sparse(nV,nV) speye(nV)/Kchrg;
                sparse(nV,nV) sparse(nV,nV)];
            Jv = J*v;
            flag = 0;
        end
        
        function [z, flag] = IDAPrecSolveFun(thisSolver, t, y, ~, ~, r, cj)
            nV = thisSolver.numVars;
            
            Kchrg = thisSolver.chargeScaleFactor;
            
            x1 = y(1:nV);
            %x2 = y(N+1:end);
            
            [~, ~, ~, dQ, dI] = thisSolver.stampFun(x1,t);
            
            P = [ dI sparse(nV, nV);
                -dQ*Kchrg   speye(nV)] + cj*[sparse(nV, nV) speye(nV)/Kchrg;
                sparse(nV, nV) sparse(nV, nV)];
            z = P\r;
            flag = 0;
        end
        
        function [R, flag] = IDAResFunEMS(thisSolver, t, y, yp)
            nV = thisSolver.numVars;
            nMCs = thisSolver.numMCs;
            nStates = thisSolver.numStates;
            
            Kchrg = thisSolver.chargeScaleFactor;
            
            v = y(1:nV);
            q = y(nV+1:2*nV);
            qp = yp(nV+1:2*nV);
            
            s = y(2*nV+1:2*nV+nStates);
            sp = yp(2*nV+1:2*nV+nStates);
            
            
            [Q, I, J] = thisSolver.stampFun(v,t,s);
            [rateVector, dRateVector, dRateVectorState, ~, ~, stateChangeMatrix] = thisSolver.ttcFun(v,t,s);
            
            [propensityVector] = thisSolver.propensityFun(rateVector, dRateVector, dRateVectorState, s);
            
            R1 = (qp/Kchrg + I - J);
            R2 = q - Q*Kchrg;
            R3 = sp - stateChangeMatrix*propensityVector;
            
            k=0;
            for i=1:1:nMCs
                m = thisSolver.circuit.MCs{i}.numStates;
                R3(k+m) = sum(s(k+1:k+m)) - sum(thisSolver.circuit.MCs{i}.stateVector);
                k=k+m;
            end
            
            R = [R1;
                R2;
                R3];
            flag = 0;
        end
        
        function [J, flag] = IDADJacFunEMS(thisSolver, t, y, ~, ~, cj)
            nV = thisSolver.numVars;
            nMCs = thisSolver.numMCs;
            nTransitions = thisSolver.numTransitions;
            nStates = thisSolver.numStates;
            
            Kchrg = thisSolver.chargeScaleFactor;
            
            v = y(1:nV);
            %q = y(nV+1:2*nV);
            s = y(2*nV+1:2*nV+nStates);
            
            [~, ~, ~, dQ, dI] = thisSolver.stampFun(v,t,s);
            
            [rateVector, dRateVector, dRateVectorState, dQtr, dItr, stateChangeMatrix] = ...
                thisSolver.ttcFun(v,t,s);
            
            [~, dPropensityVector, dPropensityVectorStatePart1, dPropensityVectorStatePart2] = ...
                thisSolver.propensityFun(rateVector, dRateVector, dRateVectorState, s);
            
            dPropensityVectorState = dPropensityVectorStatePart1 + dPropensityVectorStatePart2;
            
            J11 = dI;
            J12 = sparse(nV,nV);
            J13 = dItr;
            
            J21 = -dQ*Kchrg;
            J22 = speye(nV);
            J23 = -dQtr*Kchrg;
            
            J31 = -stateChangeMatrix*dPropensityVector;
            J32 = sparse(nStates,nV);
            J33 = -stateChangeMatrix*dPropensityVectorState;
            
            k=0;
            H = speye(nStates,nStates);
            for i=1:1:nMCs
                m = thisSolver.circuit.MCs{i}.numStates;
                J31(k+m,:) = zeros(nV,1);
                J33(k+m,:) = zeros(nStates,1);
                J33(k+m,k+1:k+m) = ones(m,1);
                H(k+m,k+m) = 0;
                k=k+m;
            end
            
            J = full([J11, J12, J13;
                J21, J22, J23;
                J31, J32, J33] + cj*[sparse(nV,nV),  speye(nV)/Kchrg,  sparse(nV,nStates);
                sparse(nV,nV),  sparse(nV,nV),    sparse(nV,nStates);
                sparse(nStates,nV), sparse(nStates,nV),   H]);
            
            flag = 0;
            j=0;
        end
        
        function [Jv, flag] = IDAJacTimesVecFunEMS(thisSolver, t, y, ~, ~, vec, cj)
            nV = thisSolver.numVars;
            nMCs = thisSolver.numMCs;
            nTransitions = thisSolver.numTransitions;
            nStates = thisSolver.numStates;
            
            Kchrg = thisSolver.chargeScaleFactor;
            
            v = y(1:nV);
            %q = y(nV+1:2*nV);
            s = y(2*nV+1:2*nV+nStates);
            
            [~, ~, ~, dQ, dI] = thisSolver.stampFun(v,t,s);
            
            [rateVector, dRateVector, dRateVectorState, dQtr, dItr, stateChangeMatrix] = ...
                thisSolver.ttcFun(v,t,s);
            
            [~, dPropensityVector, dPropensityVectorStatePart1, dPropensityVectorStatePart2] = ...
                thisSolver.propensityFun(rateVector, dRateVector, dRateVectorState, s);
            
            dPropensityVectorState = dPropensityVectorStatePart1 + dPropensityVectorStatePart2;
            
            
            J11 = dI;
            J12 = sparse(nV,nV);
            J13 = dItr;
            
            J21 = -dQ*Kchrg;
            J22 = speye(nV);
            J23 = -dQtr*Kchrg;
            
            J31 = -stateChangeMatrix*dPropensityVector;
            J32 = sparse(nStates,nV);
            J33 = -stateChangeMatrix*dPropensityVectorState;
            
            k=0;
            H = speye(nStates,nStates);
            for i=1:1:nMCs
                m = thisSolver.circuit.MCs{i}.numStates;
                J31(k+m,:) = zeros(nV,1);
                J33(k+m,:) = zeros(nStates,1);
                J33(k+m,k+1:k+m) = ones(m,1);
                H(k+m,k+m) = 0;
                k=k+m;
            end
            
            J = [J11, J12, J13;
                J21, J22, J23;
                J31, J32, J33] + cj*[sparse(nV,nV),  speye(nV)/Kchrg,  sparse(nV,nStates);
                sparse(nV,nV),  sparse(nV,nV),    sparse(nV,nStates);
                sparse(nStates,nV), sparse(nStates,nV),   H];
            
            Jv = J*vec;
            flag = 0;
        end
        
        function [z, flag] = IDAPrecSolveFunEMS(thisSolver, t, y, ~, ~, r, cj)
            nV = thisSolver.numVars;
            nMCs = thisSolver.numMCs;
            nTransitions = thisSolver.numTransitions;
            nStates = thisSolver.numStates;
            
            Kchrg = thisSolver.chargeScaleFactor;
            
            v = y(1:nV);
            %q = y(nV+1:2*nV);
            s = y(2*nV+1:2*nV+nStates);
            
            [~, ~, ~, dQ, dI] = thisSolver.stampFun(v,t,s);
            
            [rateVector, dRateVector, dRateVectorState, dQtr, dItr, stateChangeMatrix] = ...
                thisSolver.ttcFun(v,t,s);
            
            [~, dPropensityVector, dPropensityVectorStatePart1, dPropensityVectorStatePart2] = ...
                thisSolver.propensityFun(rateVector, dRateVector, dRateVectorState, s);
            
            dPropensityVectorState = dPropensityVectorStatePart1 + dPropensityVectorStatePart2;
            
            J11 = dI;
            J12 = sparse(nV,nV);
            J13 = dItr;
            
            J21 = -dQ*Kchrg;
            J22 = speye(nV);
            J23 = -dQtr*Kchrg;
            
            J31 = -stateChangeMatrix*dPropensityVector;
            J32 = sparse(nStates,nV);
            J33 = -stateChangeMatrix*dPropensityVectorState;
            
            k=0;
            H = speye(nStates,nStates);
            for i=1:1:nMCs
                m = thisSolver.circuit.MCs{i}.numStates;
                J31(k+m,:) = zeros(nV,1);
                J33(k+m,:) = zeros(nStates,1);
                J33(k+m,k+1:k+m) = ones(m,1);
                H(k+m,k+m) = 0;
                k=k+m;
            end
            
            J = [J11, J12, J13;
                J21, J22, J23;
                J31, J32, J33] + cj*[sparse(nV,nV),  speye(nV)/Kchrg,  sparse(nV,nStates);
                sparse(nV,nV),  sparse(nV,nV),    sparse(nV,nStates);
                sparse(nStates,nV), sparse(nStates,nV),   H];
            
            z = J\r;
            flag = 0;
        end
        
        
        function resFun = SEQ_TRPZResFun(thisSolver, t, Qprev, Iprev, Jprev, h, sallprev, propensityVectorPrev, propensityVectorEqPrev, dW)
            nV = thisSolver.numVars;
            nTr = thisSolver.numMCs;
            nStates = thisSolver.numStates;
            nTransitions = thisSolver.numTransitions;
            stateChangeMatrix = thisSolver.circuit.stateChangeMatrixi;
            stateChangeMatrixDiff = thisSolver.circuit.stateChangeMatrixiDiff;
            Kchrg = thisSolver.chargeScaleFactor;
            
            function R = res(x)
                % v... new voltages and possibly currents
                % s... new trap states
                v = x(1:nV);
                s = x(nV+1:nV+nStates);
                seq = x(nV+nStates+1:nV+2*nStates);
                
                sprev = sallprev(1:nStates);
                seqprev = sallprev(nStates+1:2*nStates);
                
                [Q, I, J] = thisSolver.stampFun(v,t,s);
                %[propensityVector] = thisSolver.ttcFun(v,t,s);
                %[propensityVectorEq] = thisSolver.ttcFun(v,t,seq);
                
                [rateVector, dRateVector, dRateVectorState] = ...
                    thisSolver.ttcFun(v,t,s);
                
                [propensityVector] = ...
                    thisSolver.propensityFun(rateVector, dRateVector, dRateVectorState, s);
                [propensityVectorEq] = ...
                    thisSolver.propensityFun(rateVector, dRateVector, dRateVectorState, seq);
                %[propensityVectorEqNow] = ...
                %thisSolver.propensityFun(rateVector, dRateVector, dRateVectorState, seqprev);
                
                DW = stateChangeMatrixDiff*(sqrt(propensityVectorEqPrev).*dW);
                
                R1 = Q - Qprev - h/2*(J - I + Jprev - Iprev);
                
                R2 = s - sprev - ...
                    (h/2)*(stateChangeMatrix*(propensityVector + propensityVectorPrev)) + DW;
                
                R3 = seq - seqprev - ...
                    (h/2)*(stateChangeMatrix*(propensityVectorEq + propensityVectorEqPrev));
                
                R = [Kchrg*R1; R2; R3];
            end
            
            resFun = @res;
        end
        
        function jacFun = SEQ_TRPZJacFun(thisSolver, t, h)
            nV = thisSolver.numVars;
            nStates = thisSolver.numStates;
            nTransitions = thisSolver.numTransitions;
            stateChangeMatrix = thisSolver.circuit.stateChangeMatrixi;
            Kchrg = thisSolver.chargeScaleFactor;
            
            function J = jac(x)
                % v... new voltages and possibly currents
                % s... new trap states
                v = x(1:nV);
                s = x(nV+1:nV+nStates);
                seq = x(nV+nStates+1:nV+2*nStates);
                
                [~, ~, ~, dQ, dI] = thisSolver.stampFun(v,t,s);
                %[~, dPropensityVector, dPropensityVectorState, dQst, dIst] = ...
                %    thisSolver.ttcFun(v,t,s);
                
                %[~, dPropensityVectorEq, dPropensityVectorStateEq] = ...
                %    thisSolver.ttcFun(v,t,seq);
                
                [rateVector, dRateVector, dRateVectorState, dQst, dIst] = ...
                    thisSolver.ttcFun(v,t,s);
                
                [~, dPropensityVector, dPropensityVectorStatePart1, dPropensityVectorStatePart2] = ...
                    thisSolver.propensityFun(rateVector, dRateVector, dRateVectorState, s);
                dPropensityVectorState =  dPropensityVectorStatePart1 + dPropensityVectorStatePart2;
                
                [~, dPropensityVectorEq, dPropensityVectorStateEqPart1, dPropensityVectorStateEqPart2] = ...
                    thisSolver.propensityFun(rateVector, dRateVector, dRateVectorState, seq);
                
                J11 = dQ + h/2*dI;
                J12 = dQst + h/2*dIst;
                J13 = zeros(nV, nStates);
                
                J21 = -(h/2)*(stateChangeMatrix*dPropensityVector);
                J22 = eye(nStates) - (h/2)*(stateChangeMatrix*dPropensityVectorState);
                J23 = zeros(nStates, nStates);
                
                J31 = -(h/2)*(stateChangeMatrix*dPropensityVectorEq);
                J32 = -(h/2)*(stateChangeMatrix*dPropensityVectorStateEqPart1);
                J33 = eye(nStates) - (h/2)*(stateChangeMatrix*dPropensityVectorStateEqPart2);
                
                
                J = [Kchrg*J11, Kchrg*J12, Kchrg*J13;
                    J21,       J22,       J23;
                    J31,       J32,       J33    ];
            end
            jacFun = @jac;
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %                    Non-linear solver                            %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [x, numIter] = newtSol(thisSolver, inFun, dInFun, x0, varargin)
            % NEWTSOL Newton's method for the solution of non-linear equations
            % inFun... function handle to the residual function
            % dInFun... Jacobian of the residual function
            % x0 initial guess
            % varargin can be 'Damped' for damped Newton's method
            atol = thisSolver.tolVec;
            rtol = thisSolver.reltol;
            restol = thisSolver.restol;
            damping = false;
            % is damping on?
            if (nargin > 4) && strcmp(varargin{1}, 'Damped')
                damping = true;
                imax = 4;
            end
            numIter = 0;
            convCheck = false;
            x = x0;
            inFunVal = inFun(x);
            dampInd = 1;
            numFunEval = 1;
            while convCheck == false
                dInFunVal = dInFun(x);
                try
                    lastwarn('');
                    deltax = - dInFunVal\inFunVal;
                    error(lastwarn);
                catch exception
                    rethrow(exception);
                end
                x_new = x + deltax;
                lsVecPrev = inFun(x_new);
                lsNormPrev = norm(lsVecPrev./restol);
                % if damping is on, do a line search
                if damping == true
                    for dampInd = 1:imax
                        lsVecCur = inFun(x + deltax/(2^dampInd));
                        lsNormCur = norm(lsVecCur./restol);
                        if lsNormCur > lsNormPrev
                            break;
                        end
                        lsVecPrev = lsVecCur;
                        lsNormPrev = lsNormCur;
                    end
                    alpha = 1/(2^(dampInd-1));
                    x_new = x + alpha*deltax;
                end
                numFunEval = numFunEval + dampInd;
                inFunVal = lsVecPrev;
                convCheck = (norm((x_new - x)./(atol + rtol*abs(x_new))) <= 1) && ...
                    (lsNormPrev <= 1);
                x = x_new;
                numIter = numIter + 1;
            end
        end
        
        
        function [propensityVector, dPropensityVector, dPropensityVectorStatePart1, dPropensityVectorStatePart2] ...
                = propensityFun(thisSolver, rateVector, dRateVector, dRateVectorState, s)
            
            m=0;
            n=0;
            k=m;
            l=n;
            thisCircuitMC = thisSolver.circuit;
            numCircuitStates = thisCircuitMC.numStates;
            nMCs = thisCircuitMC.numMCs;
            parentComponentCurrent = thisCircuitMC.MCs{1}.parentComponent;
            statesCircComp = parentComponentCurrent.numStates;
            transitionsCircComp = parentComponentCurrent.numTransitions;
            for j=1:nMCs
                
                MC = thisCircuitMC.MCs{j};
                
                if isequal(parentComponentCurrent, MC.parentComponent) == false
                    parentComponentPrev = parentComponentCurrent;
                    parentComponentCurrent = MC.parentComponent;
                    m = m + parentComponentPrev.numTransitions;
                    n = n + parentComponentPrev.numStates;
                    k = m;
                    l = n;
                end
                
                statesMC = MC.numStates;
                transitionsMC = MC.numTransitions;
                
                [sPropensityVector, sdPropensityVector, sdPropensityVectorStatePart1, sdPropensityVectorStatePart2] ...
                    = MC.propensityEval(rateVector(k+1:k+transitionsMC,1), dRateVector(k+1:k+transitionsMC,:), dRateVectorState(k+1:k+transitionsMC,:), s(l+1:l+statesMC,1), l, numCircuitStates);
                
                propensityVector(k+1:k+transitionsMC,1) = sPropensityVector;
                dPropensityVector(k+1:k+transitionsMC,:) = sdPropensityVector;
                dPropensityVectorStatePart1(k+1:k+transitionsMC,:) = sdPropensityVectorStatePart1;
                dPropensityVectorStatePart2(k+1:k+transitionsMC,:) = sdPropensityVectorStatePart2;
                k=k+transitionsMC;
                l=l+statesMC;
            end
        end
        
        %end methods
    end
end