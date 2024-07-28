classdef lptvMC < handle
	% This class takes a harmonic balance solver object as input and
	% calculates steady-state timing jitter variance slope as output.
    
    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------
 
	properties (SetAccess = private)
		circuit = circuit.empty;
		hbSolver = harmonicBalanceSolverMC.empty;
		Tf = 0;			% fundamental period;
		nTimeStep = 0;	% number of time steps
		tVec;
		nVar = 0;
        nMCs = 0;
        nStates = 0;
        nTransitions = 0;
		CF;		
		GF;		
		BF; 
		Yss;
		JM;
		dQM;
		dIM;
		dBM;
		LM;
		DN;
		DNinv;
		OM;
		yIdx;	% index of the circuit variable of interest
		autonom = false;
		Un; % nullspace for GF
		Vn; % nullspace for GF'
		jitterSlope = 0;
		genEigU = []; % solution of the generalized eigenvalue problem for the lptv system
		genEigV = [];
		genEigDu = [];
		genEigDv = [];
	end

	methods (Access = private)

	end

	methods

		function lptvGen(thisLptv, slvr)

			if isa(slvr, 'harmonicBalanceSolverMC') &&...
								   ~isempty(slvr.Y) &&...
								   ~isempty(slvr.JM) &&...
								   ~isempty(slvr.dQM) &&...
								   ~isempty(slvr.dIM) &&...
								   ~isempty(slvr.dBM)
				
				thisLptv.hbSolver = slvr;
				thisLptv.nTimeStep = slvr.nFreq*2+1;
				thisLptv.nVar = slvr.nVar + slvr.nStates - slvr.nMCs;
				thisLptv.nMCs = slvr.nMCs;
                thisLptv.nStates = slvr.nStates;
                thisLptv.nTransitions = slvr.nTransitions;

				thisLptv.Yss = slvr.Y;
				thisLptv.JM = slvr.JM;
				thisLptv.dQM = slvr.dQM;
				thisLptv.dIM = slvr.dIM;
				thisLptv.dBM = slvr.dBM;
				thisLptv.DN = slvr.DN;
				thisLptv.DNinv = slvr.DNinv;
				thisLptv.OM = slvr.OM;
				thisLptv.tVec = slvr.tVec;
				thisLptv.Tf = slvr.Tf;

				thisLptv.autonom = slvr.autonom;
			else
				error('One or more properties of the harmonicBalanceSolver are empty');
            end
            
			M = thisLptv.nTimeStep;
			P = thisLptv.nTransitions;

			[~, DPinv] = dftmtxNM(P,M);

			OM = thisLptv.OM;
			DN = thisLptv.DN;
			dQM = blkdiag(thisLptv.dQM{:});
			dIM = blkdiag(thisLptv.dIM{:});
			DNinv = thisLptv.DNinv;
			dBM = blkdiag(thisLptv.dBM{:});


			% 258.54 secs
			CF = DN*dQM*DNinv;
			GF = OM*CF + DN*dIM*DNinv;

			thisLptv.GF = GF;
			thisLptv.CF = CF;
            
			if thisLptv.nMCs > 0
				BF = DN*dBM*DPinv;
				thisLptv.BF = BF;
			end
		end

		function [un, vn] = floquetVector(thisLptv,yIdx)
            %keyboard
			N = thisLptv.nVar;
			K = (thisLptv.nTimeStep-1)/2;
			un = thisLptv.OM*(thisLptv.hbSolver.YF);
			% TODO: implement a check here for the 'singular' eigenvalue
			[vn, ~] = eigs(thisLptv.GF', 1, 'sm'); 
			% make the elements of vn conjugate symmetric so that its inverse Fourier
			% transform is real
			vn = vn * exp(-1i*angle(vn(N*(K-1)+yIdx)/conj(vn(N*(K+1)+yIdx)))/2);
			% scale vn such that its bi-CF orthonormal with un
			vn = vn/(vn'*thisLptv.CF*un);
			%vn = vn/(vn'*un);
			thisLptv.Un = un;
			thisLptv.Vn = vn;
		end

		function js = computeJS(thisLptv, yIdx)
			N = thisLptv.nVar;
			M = thisLptv.nTimeStep;
            T = thisLptv.nTransitions;
            if thisLptv.jitterSlope ==0 
				if isempty(thisLptv.Vn) || isempty(thisLptv.Un)
					[~, vn] = thisLptv.floquetVector(yIdx);
				else
					vn = thisLptv.Vn;
                end
                
				vnt = reshape(real(ifftNM(vn, N, M)), N, M);
				vt = zeros(T,M);
				vtTvt = zeros(1,M);
                
				for i=1:thisLptv.nTimeStep
					vt(:,i) = (vnt(:,i)'*thisLptv.dBM{i})';
					vtTvt(i) = vt(:,i)'*vt(:,i);
				end

				js = trapz(thisLptv.tVec, vtTvt)/thisLptv.Tf;
				thisLptv.jitterSlope = js;
			else
				js = thisLptv.jitterSlope;
			end
        end

	% end methods
	end
end
