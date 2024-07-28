classdef harmonicBalanceSolverMC < handle
	% HARMONICBALANCESOLVERMC harmonic balance periodic steady-state analysis.
	% Works for autonomous neuronal circuits. A fundamental period estimate 
    % Tf_est should be provided.
    
    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------

	properties (Constant)
        keyNames = {'currentVariables'
					'transientSolverMC'
					'yss_est'
					'Tf_est'		% fundamental period estimate
					'Tf'
					'nFreq'
					't0'
                    'reltol'
                    'abstol_v'
                    'abstol_c'
					'yIdx'			% index of he periodic variable in the circuit
					'v0'			% start value for the waveform
                    'linSolverType'
					'gmrestart'
					'gpufft'
                    };

        validOnOff = {'on', 'off'};
        validLinSolverTypes = {'BSLASH', 'GMRES'};

	end

	properties (SetAccess = private)
		transientSolverMC = transientSolverMC.empty;
		circuit = circuit.empty;
		yss_est;		% estimate of steady state solution
		Tf_est = 0;			% fundamental period estimate for autonomous systems  
		yss;
		Tf = 0;				% fundamental period -- IF NOT GIVEN THE CIRCUIT IS ASSUMED AUTONOMOUS
		fnconst = 0;		% frequency normalizing constant.
		t0 = inf;
		tVec;
		nFreq = 0;			% number of frequencies -- same as number of time steps
		fVec;
		currentVariables = 'off';
		stampFun;
		ttcFun;
		nVar = 0;			% total number of variables 
		nVoltVar = 0;
		nCurrVar = 0;
        nMCs = 0;
        nStates = 0;
        nTransitions = 0;
        stateChangeMatrix = 0;
		reltol;
		abstol_v;
		abstol_c;
		tolVec;
		OM;					% frequency domain differentiation matrix diagonal NMxNM
		DN;					% dFT matrix NMxNM
		DNinv;				% inverse DFT matrix NMxNM
		JF;					% fourier transform of NMx1 source vector in form of [J(t_1);...;J(t_M)]
		Y;					% solution array
		YF;					% solution in Fourier domain
		dQM;
		dIM;
		dBM;
		JM;
		autonom = false;	% is the circuit an autonomous oscillator
		resJacFun;			% residual and jacobian function -- changes acc. to autonom or not
		resFun;
		jacFun;
		yIdx = 0;			% index of the node of interest
		v0 = inf;			% start value for autonomous circuits, to fix the phase
		linSolverType = 'GMRES';
		gmresFlag = false;
		gmrestart = 100;
		gpufft = false;
	end

	methods
		function obj = harmonicBalanceSolverMC(ckt, varargin)
			if nargin > 0
                if nargin > 1
                    obj.setOptions(varargin{:});
                end
				obj.setCircuit(ckt);
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
            
            if all(strcmp(thisSolver.currentVariables, harmonicBalanceSolverMC.validOnOff) == 0)
                msgStr = sprintf('Unrecognized current variables option: %s', thisSolver.currentVariables);
                msgStr = sprintf('%s\nValid choices are:', msgStr);
                for i=1:numel(harmonicBalanceSolverMC.validOnOff)
                    msgStr = sprintf('%s\n%d. %s', msgStr, i, harmonicBalanceSolverMC.validOnOff{i});
                end
                error(msgStr);
            end

            if all(strcmp(thisSolver.linSolverType, harmonicBalanceSolverMC.validLinSolverTypes) == 0)
                msgStr = sprintf('Unrecognized linear solver type: %s', thisSolver.linSolverType);
                msgStr = sprintf('%s\nValid choices are:', msgStr);
                for i=1:numel(harmonicBalanceSolverMC.validLinSolverTypes)
                    msgStr = sprintf('%s\n%d. %s', msgStr, i, harmonicBalanceSolverMC.validLinSolverTypes{i});
                end
                error(msgStr);
            end

			if strcmp(thisSolver.linSolverType, 'GMRES')
				thisSolver.gmresFlag = true;
			end
			
			if thisSolver.Tf == 0
				thisSolver.autonom = true;
				if thisSolver.gmresFlag == true
					thisSolver.resJacFun = @thisSolver.resJacTimesVecFunAutonom;
				else
					thisSolver.resJacFun = @thisSolver.resJacFunAutonom;
				end
			else
				thisSolver.autonom = false;
				thisSolver.Tf_est = thisSolver.Tf;
				thisSolver.resJacFun = @thisSolver.resJacFunForced;
				if thisSolver.gmresFlag == true
					thisSolver.resJacFun = @thisSolver.resJacTimesVecFunForced;
				else
					thisSolver.resJacFun = @thisSolver.resJacFunForced;
				end
			end

			if ~isempty(thisSolver.yss_est) && isinf(thisSolver.t0)
				error('You provided a steady state estimate but not a start time!')
			end

			if ~isempty(thisSolver.yss_est) && thisSolver.Tf_est == 0
				error('You provided a steady state estimate but not a fundamental period!')
			end

			if ~isempty(thisSolver.tVec) &&...
			    thisSolver.nFreq ~= 0 &&...
			    (2*thisSolver.nFreq + 1) ~= length(thisSolver.tVec)
				error('Size of tVec must match nFreq, length(tVec) = 2*nFreq + 1!');
			end

           
        end

		function setCircuit(thisSolver, ckt)
			thisSolver.circuit = ckt;
			thisSolver.nMCs = ckt.numMCs;
            thisSolver.nStates = ckt.numStates;
            thisSolver.nTransitions = ckt.numTransitions;
            thisSolver.stateChangeMatrix = ckt.stateChangeMatrixi;
			nMCs = thisSolver.nMCs;
            nStates = thisSolver.nStates;

            switch thisSolver.currentVariables
                case 'on'
                    thisSolver.nVar = ckt.numVars;
                    thisSolver.nVoltVar = ckt.numIndepVoltVars;
					thisSolver.nCurrVar = ckt.numCurrentVars;
					thisSolver.stampFun = @ckt.stamp;
					if nMCs ~= 0
						thisSolver.ttcFun = @ckt.ttc;
					end
                case 'off'
                    thisSolver.nVar = ckt.numVarsnc;
                    thisSolver.nVoltVar = ckt.numVarsnc;
					thisSolver.nCurrVar = 0;
					thisSolver.stampFun = @ckt.stampnc;
					if nMCs ~= 0
						thisSolver.ttcFun = @ckt.ttcnc;
					end
            end         

            % construct the tolerance vector
            nV = thisSolver.nVar;
            nVV = thisSolver.nVoltVar;
			nCV = thisSolver.nCurrVar;
            
            tolVec = [ones(nVV,1).*thisSolver.abstol_v;
					  ones(nCV,1).*thisSolver.abstol_c;
					  ones(nStates-nMCs,1).*thisSolver.abstol_v];

			% replicate this tolerance vector M times
			thisSolver.tolVec = tolVec((1:nV+nStates-nMCs)'*ones(1,thisSolver.nFreq*2+1),:);

			% if the circuit is autonomous, then we need an additional equation
			if thisSolver.autonom == true
				thisSolver.tolVec = [thisSolver.tolVec; thisSolver.abstol_v];
			end

                    
			if ~isempty(thisSolver.yss_est) && ...
				size(thisSolver.yss_est, 1) ~= thisSolver.nVar + thisSolver.nStates
				error('Number of variables in yss_est and the circuit are different!')
			end

		end

		function setPeriod(thisSolver, T)
			thisSolver.Tf = T;
		end

		function solve(thisSolver)
			N = thisSolver.nVar + thisSolver.nStates - thisSolver.nMCs;
			M = thisSolver.nFreq*2 + 1;

			if isempty(thisSolver.yss_est)
				error('Please run estimateSS method first or provide a steady state solution estimate');
			end

			thisSolver.operatorGen();

			x0 = thisSolver.yss_est(:);
			
			XF0 = fftNM(x0,N,M,thisSolver.gpufft);
			if thisSolver.autonom == true
				XF0 = [XF0; 1/thisSolver.Tf_est/thisSolver.fnconst];
			end

			disp('Starting harmonic balance solution');
            
			[YF, ~] = thisSolver.newtSol(thisSolver.resJacFun, XF0);
			disp('Convergence.');

			if thisSolver.autonom == true
                %keyboard
				%thisSolver.Tf = real(1/thisSolver.fnconst/YF(N*M + 1));
                thisSolver.Tf = 1/real(YF(N*M+1)*thisSolver.fnconst);
				YF = YF(1:N*M);
            end
            
			% NOTE: here we do not explicitely save the Jacobians etc.
			% for the lptv system. In the newton solver, these variables are saved at each step.
			% so for this solver, there is no reason to repeat a equation evaluation routine here.
			thisSolver.YF = YF;
			thisSolver.Y = reshape(real(ifftNM(YF,N,M,thisSolver.gpufft)), N, thisSolver.nFreq*2+1);
		end

		function estimateSS(thisSolver, tSlvr)
			% estimate steady state periodic solution
			% tSlvr...transient solver object
			% yIdx...index of the variable of interest
			% There is a logic here that acts differently according to 
			% whether t0 and/or Tf is provided

			if strcmp(thisSolver.currentVariables, tSlvr.currentVariables) == 0
				error(['The currentVariables option of the transient solver '...
					   'does not match with this option in the harmonicBalanceSolver']);
			end

			if thisSolver.circuit ~= tSlvr.circuit
				error('Please provide a transientSolver for the same circuit.');
			end


			if isempty(tSlvr.Y) || isempty(tSlvr.T)
				disp('Starting transient simulation.')
				tSlvr.solve();
			else
				disp('Skipping transient simulation because the provided solver already contains a solution');
            end
            
            thisSolver.transientSolverMC = tSlvr;
			t_slvr = tSlvr.T;
			y_est = tSlvr.Y(thisSolver.yIdx,:);

			[t_beg, t_end] = thisSolver.findTBegEnd(t_slvr, y_est);
			thisSolver.t0 = t_beg;
			Tf = t_end-t_beg;

			if thisSolver.Tf_est == 0
				thisSolver.Tf_est = Tf;
			end

			thisSolver.fnconst = 1/thisSolver.Tf_est/1e2;

			% create time vector 
			thisSolver.tVecGen();
            
            %exclude one state variable for each MC
            nMCs = thisSolver.nMCs;
            excludedState = 1;
            excludeStateInd(1) = excludedState;
            for ind=2:1:nMCs
                m = thisSolver.circuit.MCs{ind-1}.numStates;
                excludeStateInd(ind) = excludeStateInd(ind-1) + m;
            end
            S = tSlvr.S;
            S(excludeStateInd,:)=[];

			nV = tSlvr.numVars;
			thisSolver.yss_est = interp1(t_slvr.',[tSlvr.Y(1:nV,:); S].', thisSolver.tVec.').';

			if thisSolver.autonom == true && isinf(thisSolver.v0)
				thisSolver.v0 = thisSolver.yss_est(thisSolver.yIdx,1);
			end
		end

		function [t_beg, t_end] = findTBegEnd(thisSolver, t_slvr, y_est)
			% find start and end times for the periodic ss estimate
			% for four different cases depending on whether (t0, Tf_est)
			% are provided

			if ~isinf(thisSolver.t0)
				t_beg = thisSolver.t0;
				v_beg = interp1(t_slvr,y_est, t_beg);

				if isnan(v_beg)
					error('the transient simulation does not include t0!');
				end

				if thisSolver.Tf_est ~= 0;
					t_end = t_beg + thisSolver.Tf_est;

					if t_slvr(end) < t_end
						error('The fundamental period you have provided seems to be too long.');
					end
				else
					y_estShift = y_est - v_beg;
					% calculate crossing points. We have always: y_est(xIdx) <= v_beg 
					xPointVec = sign(y_estShift(1:end-1)) .* sign(y_estShift(2:end)); 
					xIdx = find(xPointVec ~= 1);

					if any(diff(xIdx) < 5)
						error(['Points with f(t) = f(t0) are found ',...
								 'that are close to each other. ',...
								 'Maybe you shoud specify another t0']);
					end

					% find the index of the beginning time
					[~, begIdx] = min(abs(t_slvr-t_beg));
					% always extract the time index before the crossing
					begIdx = begIdx - 1*(prod(y_estShift(begIdx:begIdx+1)) > 0);

					if find(begIdx == xIdx) > length(xIdx) - 2;
						error(['Could not find a third point with f(t) = f(t0). ',...
							   'Maybe the transient simulation does not ',...
							   'continue for one period after t0.']);
					end

					endIdx = xIdx(find(xIdx == begIdx) + 2);
					t_end = t_slvr(endIdx) + y_estShift(endIdx)/...
							(y_estShift(endIdx)-y_estShift(endIdx+1)) * ...
							(t_slvr(endIdx+1) - t_slvr(endIdx));
				end
			else
				midVal = (max(y_est) + min(y_est))/2;
				y_estShift = y_est - midVal;
				
				% crossing points
				xPointVec = sign(y_estShift(1:end-1)) .* sign(y_estShift(2:end)); 
				xIdx = find(xPointVec ~= 1);

				if any(diff(xIdx) < 5)
					warning(['Points with f(t) = f(t0) are found ',...
							 'that are close to each other. ',...
							 'Maybe you shoud specify t0']);
				end
				
				% linear interpolation
				begIdx = xIdx(end-2);
				t_beg = t_slvr(begIdx) + y_estShift(begIdx)/...
						(y_estShift(begIdx)-y_estShift(begIdx+1)) * ...
						(t_slvr(begIdx+1) - t_slvr(begIdx));

				if thisSolver.Tf_est ~= 0
					t_end = t_beg + thisSolver.Tf_est;
					if t_slvr(end) < t_end
						error('The fundamental period you have provided seems to be too long.');
					end
				else 
					endIdx = xIdx(end);
					t_end = t_slvr(endIdx) + y_estShift(endIdx)/...
							(y_estShift(endIdx)-y_estShift(endIdx+1)) * ...
							(t_slvr(endIdx+1) - t_slvr(endIdx));
				end
			end
		end

		function operatorGen(thisSolver)
			thisSolver.tVecGen();
			thisSolver.DNGen();
			thisSolver.OMGen();
			thisSolver.JFGen();
		end

		function tVecGen(thisSolver)
			K = thisSolver.nFreq; % number of freqs
			M = 2*K + 1;			% number of data points
			Tf = thisSolver.Tf_est;
			t0 = thisSolver.t0;

			% the last point is the same with the first; ignore it.
			thisSolver.tVec = cumsum([t0, (Tf/M)*ones(1,M-1)]);
		end

		function DNGen(thisSolver)
			K = thisSolver.nFreq; % number of freqs
			M = 2*K + 1;			% number of data points
			N = thisSolver.nVar + thisSolver.nStates - thisSolver.nMCs;	% number of variables in the circuit eqs
			Tf = thisSolver.Tf_est;
			fc = 1/Tf;
			%tVec = thisSolver.tVec;

			thisSolver.fVec = fc*(-K:K)/M;

			[thisSolver.DN, thisSolver.DNinv] = dftmtxNM(N,M);
		end

		function OMGen(thisSolver)
			K = thisSolver.nFreq; % number of freqs
			M = 2*K + 1;			% number of data points
			N = thisSolver.nVar + thisSolver.nStates - thisSolver.nMCs;	% number of variables in the circuit eqs
			Tf = thisSolver.Tf_est;
			fc = 1/Tf;
			
			% generate OM -- freq domain differentiation matrix
			diagVec = 1i*2*pi*fc*(-K:K);
			diagVec = diagVec(ones(N,1),:);
			OM = sparse(1:N*M,1:N*M,diagVec);

			thisSolver.OM = OM;
		end

		function JFGen(thisSolver)
			nEV = thisSolver.nVar;
			nMCs = thisSolver.nMCs;
            nStates = thisSolver.nStates;
			tVec = thisSolver.tVec;

			K = thisSolver.nFreq; % number of freqs
			M = 2*K + 1;			% number of data points
			N = nEV + nStates - nMCs;	% number of variables in the circuit eqs
			
			% generate JF -- source vector
			J = zeros(N,M);
			sFun = thisSolver.stampFun;
			parfor i=1:M
				[~, ~, Ji] = sFun(zeros(nEV,1), tVec(i));
				% NOTE: room for improvement -- do not compute FT of 
				% the zeros below inserted for the MC variables.
				J(:,i) = [Ji; zeros(nStates-nMCs,1)];
			end

			thisSolver.JM = J;
			thisSolver.JF = fftNM(J(:),N,M,thisSolver.gpufft);
		end

		function [QM, IM, dQM, dIM] = stampNM(thisSolver, x)
			% x is NxM here
			sFun = thisSolver.stampFun;
			tFun = thisSolver.ttcFun;
			N = thisSolver.nVar + thisSolver.nStates - thisSolver.nMCs;
			K = thisSolver.nFreq;
			M = 2*K+1;
			nMCs = thisSolver.nMCs;
            nStates = thisSolver.nStates;
            nTransitions = thisSolver.nTransitions;
            stateChangeMatrix = thisSolver.stateChangeMatrix;
			nEV = thisSolver.nVar; % number of electrical variables
			tVec = thisSolver.tVec;
			snconst = 1e10/thisSolver.Tf_est;
            
            excludedState = 1;
            excludeStateInd(1) = excludedState;
            for ind=2:1:nMCs
                m = thisSolver.circuit.MCs{ind-1}.numStates;
                excludeStateInd(ind) = excludeStateInd(ind-1) + m;
            end

			QM = zeros(N, M);
			IM = zeros(N, M);
			dQM = cell(1,M);
			dIM = cell(1,M);
			dBM = cell(1, M);

			v = x(1:nEV,:);
			s = x(nEV+1:N,:);
			% XFRS = reshpae(XF,N,M);
			% S = XFRS(N-nTr+1:N,:);
			% S = S(:);

			% NOTE: room for improvement -- you are computing the inverse FT
			% of the MC variables and then again the FT below.
			% For now I am leaving this as is because I don't want to change the
			% order of variables, i.e. separate the electrical and MC variables
			% in the big XF vector. The same applies for the Jacobians.
			parfor i=1:M
				% NOTE: the external forcing, J, is not computed here.
				vi = v(:,i);
				if nMCs == 0
					[Qi, Ii, ~, dQi, dIi] = sFun(vi, tVec(i));

					QM(:,i) = Qi;
					IM(:,i) = Ii;
					dQM{i} = dQi;
					dIM{i} = dIi;
				else
					si = s(:,i);
                    
                    %reconstruct state vector si with excluded states
                    k=0;
                    for j=1:1:nMCs
                        insertInd = excludeStateInd(j);
                        m = thisSolver.circuit.MCs{j}.numStates;
                        foundState = sum(thisSolver.circuit.MCs{j}.stateVector) - sum(si(k+1:k+m-1));
                        si = [si(1:insertInd-1); foundState; si(insertInd:end)];
                        k=k+m;
                    end
                    
                    [Qi, Ii, ~, dQi, dIi] = sFun(vi, tVec(i), si);
                    [rateVector, dRateVector, dRateVectorState, dQst, dIst] = tFun(vi, tVec(i), si);  
                    [propensityVector, dPropensityVector, dPropensityVectorStatePart1, dPropensityVectorStatePart2] = ...
                        thisSolver.pFun(rateVector, dRateVector, dRateVectorState, si);
                    dPropensityVectorState =  dPropensityVectorStatePart1 + dPropensityVectorStatePart2;

                    QM2 = si/snconst;
                    IM2 = -stateChangeMatrix*propensityVector/snconst;
                    dBM2 = stateChangeMatrix.*repmat(sqrt(propensityVector)',nStates,1)/snconst;
                    dQM22 = speye(nStates)/snconst;
                    dIM21 = -stateChangeMatrix*dPropensityVector/snconst;
                    dIM22 = -stateChangeMatrix*dPropensityVectorState/snconst;
                                        
                    QM_ = [Qi; QM2];
                    IM_ = [Ii; IM2];
                    dQM_ = [dQi, dQst;
                              zeros(nStates,nEV)/snconst, dQM22];
                    dIM_ = [dIi, dIst;
                              dIM21, dIM22];
                    dBM_ = [zeros(nEV, nTransitions); dBM2];
                    
                    %exclude one state variable for each MC
                    stateInd=0;
                    for j=1:nMCs
                        m = thisSolver.circuit.MCs{j}.numStates;
                        excludeColumndQ = -dQM_(:,excludeStateInd(j)+nEV);
                        excludeColumndI = -dIM_(:,excludeStateInd(j)+nEV);
                        insertMatrixdQ = repmat(excludeColumndQ,1,m);
                        insertMatrixdI = repmat(excludeColumndI,1,m);
                        dQM_(:,nEV+stateInd+1:nEV+stateInd+m) = dQM_(:,nEV+stateInd+1:nEV+stateInd+m) + insertMatrixdQ;
                        dIM_(:,nEV+stateInd+1:nEV+stateInd+m) = dIM_(:,nEV+stateInd+1:nEV+stateInd+m) + insertMatrixdI;
                        stateInd = stateInd + m;
                    end
                    QM_(excludeStateInd+nEV,:)=[];
                    IM_(excludeStateInd+nEV,:)=[];
                    dQM_(:,excludeStateInd+nEV)=[];
                    dQM_(excludeStateInd+nEV,:)=[];
                    dIM_(:,excludeStateInd+nEV)=[];
                    dIM_(excludeStateInd+nEV,:)=[];
                    dBM_(excludeStateInd+nEV,:)=[];
                    
                    %save reduced vectors and matrices
                    QM(:,i) = QM_;
                    IM(:,i) = IM_;
                    dQM{i} = dQM_;
                    dIM{i} = dIM_;
                    dBM{i} = sparse(dBM_);
				end
			end
			% ATTENTION: the properties below are saved in order to be used with
			% an lptv. Do not use them as global variables in the class internal calculations.
			% Use the local output variables for that purpose.
			thisSolver.dQM = dQM;
			thisSolver.dIM = dIM;
			thisSolver.dBM = dBM;
		% end stampNC
		end

		function [R, dR, x, R1] = resJacFunForced(thisSolver, XF)
			DN = thisSolver.DN;
			OM = thisSolver.OM;
			DNinv = thisSolver.DNinv;
			JF = thisSolver.JF;
			N = thisSolver.nVar + thisSolver.nStates - thisSolver.nMCs;
			K = thisSolver.nFreq;
			M = 2*K+1;

			x = reshape(real(ifftNM(XF,N,M,thisSolver.gpufft)), N, M);

			[QM, IM, dQM, dIM] = thisSolver.stampNM(x);

			QM = QM(:);
			IM = IM(:);
			dQM = blkdiag(dQM{:});
			dIM = blkdiag(dIM{:});

			R1 = OM*DN*QM;
			R2 = DN*IM;
			R = R1 + R2 - JF;
			dR = OM*DN*dQM*DNinv + DN*dIM*DNinv;
		end

		function [R, dR, dummyOut] = resJacFunAutonom(thisSolver, XF)
			N = thisSolver.nVar + thisSolver.nStates - thisSolver.nMCs;
			K = thisSolver.nFreq;
			M = 2*K+1;
			NM = N*M;
			yIdx = thisSolver.yIdx;
			fnconst = thisSolver.fnconst;
			v0 = thisSolver.v0;

			f_est = real(XF(NM+1)*fnconst);
			Tf_est = 1/f_est;
			thisSolver.Tf_est = Tf_est;

			thisSolver.tVecGen();
			thisSolver.OMGen();

			[R, dR, x, R1] = thisSolver.resJacFunForced(XF(1:NM));

			R = [R; x(yIdx,1)-v0];
			%R = [R; imag(XF((K+1)*N+yIdx))];

			dR22 = 0;

			dR12 = R1/f_est*fnconst;
			dR21 = sparse(1, yIdx:N:NM, thisSolver.DNinv(yIdx,yIdx:N:NM), 1, NM); % all ones
			%dR21 = sparse(1, (K+1)*N+yIdx, 1i/2, 1, NM);
			dR = [dR,	dR12;
				  dR21, dR22];
			dummyOut = [];
		end

		function [R, JxVFun, precSolveFun, x, R1] = resJacTimesVecFunForced(thisSolver, XF)
			OM = thisSolver.OM;
			JF = thisSolver.JF;
			N = thisSolver.nVar + thisSolver.nStates - thisSolver.nMCs;
			K = thisSolver.nFreq;
			M = 2*K+1;
			NM = N*M;
			j2piKfVec = OM(1:(N*M+1)*N:NM^2);


			x = reshape(real(ifftNM(XF,N,M,thisSolver.gpufft)), N, M);
			[QM, IM, dQM, dIM] = thisSolver.stampNM(x);

			dQMavg = sparse(N,N);
			dIMavg = sparse(N,N);
			for i=1:M
				dQMavg = dQMavg + dQM{i};
				dIMavg = dIMavg + dIM{i};
            end
    
			QM = QM(:);
			IM = IM(:);
			dQM = blkdiag(dQM{:});
			dIM = blkdiag(dIM{:});

			R1 = OM*fftNM(QM,N,M,thisSolver.gpufft);
			R = R1 + fftNM(IM,N,M,thisSolver.gpufft) - JF;

			function jxv = JxV(VF)
				v = ifftNM(VF,N,M,thisSolver.gpufft);
				jxv = OM*fftNM(dQM*v,N,M,thisSolver.gpufft) +...
						fftNM(dIM*v,N,M,thisSolver.gpufft);
			end

			function X = precSolve(Y)
				Y = reshape(Y,N,M);
				X = zeros(N,M);
				parfor j=1:M
					X(:,j) = (j2piKfVec(j)*dQMavg + dIMavg)\Y(:,j);
				end
				X = X(:);
			end

			JxVFun = @JxV;
			precSolveFun = @precSolve;
		end

		function [R, JxVFun, precSolveFun] = resJacTimesVecFunAutonom(thisSolver, XF)
			N = thisSolver.nVar + thisSolver.nStates - thisSolver.nMCs;
			K = thisSolver.nFreq;
			M = 2*K+1;
			NM = N*M;
			yIdx = thisSolver.yIdx;
			fnconst = thisSolver.fnconst;
			v0 = thisSolver.v0;

			f_est = real(XF(NM+1)*fnconst);
			Tf_est = 1/f_est;
			thisSolver.Tf_est = Tf_est;

			thisSolver.tVecGen();
			thisSolver.OMGen();

			[R, JxVFunFrcd, precSolveFunFrcd, x, R1] = thisSolver.resJacTimesVecFunForced(XF(1:NM));

			R = [R; x(yIdx,1)-v0];

			dR22 = 0;

			dR12 = R1/f_est*fnconst;
			dR21 = sparse(1, yIdx:N:NM, thisSolver.DNinv(yIdx,yIdx:N:NM), 1, NM); % all ones

			function jxv = JxV(VF)
				XF = VF(1:NM);
				fc = real(VF(NM+1));
				jxv = [JxVFunFrcd(XF);
						dR21*XF]		+	 [dR12;
											  dR22]*fc;
			end

			function X = precSolve(Y)
				X = [precSolveFunFrcd(Y(1:NM)); 
						Y(NM+1)];
			end

			JxVFun = @JxV;
			precSolveFun = @precSolve;
		end

		function displaySolutionF(thisSolver, YF)
			N = thisSolver.nVar + thisSolver.nStates - thisSolver.nMCs;
			K = thisSolver.nFreq;
			M = 2*K+1;
			figure;
            plot(thisSolver.tVec, real(reshape(ifftNM(YF(1:N*M),N,M,thisSolver.gpufft),N,M)));
            %drawnow
            %pause(0.05);
		end

		function displaySolution(thisSolver)
			plot(thisSolver.tVec,thisSolver.Y);
		end

		function lptvGen(thisSolver)

		end

        function [x, numIter] = newtSol(thisSolver, resJacFun, x0)
            % NEWTSOL Newton's method for the solution of non-linear equations
			% resJacFun ... function handle to compute the residual and its Jacobian
            % x0 initial guess
            
            atol = thisSolver.tolVec;
            rtol = thisSolver.reltol;
            restol = thisSolver.abstol_c;
			gmresFlag = thisSolver.gmresFlag;
			gmrestart = thisSolver.gmrestart;
            damping = true;
            imax = 4;
                        
            numIter = 0;
            convCheck = false;
            x = x0;
            [funVal, dFun, precMtx] = resJacFun(x);
			%thisSolver.displaySolutionF(x0);

			if gmresFlag == true
				fprintf('Starting Newton''s Method with GMRES\n');
			else
				fprintf('Newton iterations: 00\n');
            end
            
            numFunEval = 1;
            while convCheck == false
                try
                    lastwarn('');
					if thisSolver.gmresFlag == true
						deltax = gmres(dFun, funVal, gmrestart, restol, [], precMtx);
                        %keyboard
                    else
						deltax = dFun\funVal;
                        %keyboard
					end
                    %error(lastwarn);
                catch exception
                    rethrow(exception);
                end
                fprintf('norm of delta_x = %e\n',norm(deltax));
                x_new = x - deltax;
				%thisSolver.displaySolutionF(x_new); 
                lsVecPrev = resJacFun(x_new);
                lsNormPrev = abs(norm(lsVecPrev./restol));
                % if damping is on, do a line search
                if damping == true
                    for dampInd = 1:imax
                        lsVecCur = resJacFun(x - deltax/(2^dampInd));
                        lsNormCur = abs(norm(lsVecCur./restol));
                        if lsNormCur > lsNormPrev
                            break;
                        end
                        lsVecPrev = lsVecCur;
                        lsNormPrev = lsNormCur;
                    end
                    alpha = 1/(2^(dampInd-1));
                    fprintf('alpha = %e\n',alpha);
                    x_new = x - alpha*deltax;
                end
                %thisSolver.displaySolutionF(x_new);
                numFunEval = numFunEval + imax;
                
                [funVal, dFun, precMtx] = resJacFun(x_new);
                
                convCheck = (norm((x_new - x)./(atol + rtol*abs(x_new))) <= 1) && ...
                            (norm(funVal)/restol <= 1);
				fprintf('step size = %e\n',(norm((x_new - x)./(atol + rtol*abs(x_new)))));
				fprintf('function value = %e\n',(norm(funVal)/restol));
                x = x_new;
                numIter = numIter + 1;

				if gmresFlag == false
					fprintf('Newton iterations: %.2i\n\n', numIter);
				end
            end
            
        end
        
        function [propensityVector, dPropensityVector, dPropensityVectorStatePart1, dPropensityVectorStatePart2] ...
                = pFun(thisSolver, rateVector, dRateVector, dRateVectorState, s)
            
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
	% end methods
	end
% end classdef
end
