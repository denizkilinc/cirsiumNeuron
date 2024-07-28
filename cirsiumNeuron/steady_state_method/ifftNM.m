function x = ifftNM(XF, N, M, gpuFlag)
    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------
    
	% N... number of variables
	% M... number of time points
	K = (M-1)/2; % number of freqs
	nCol = size(XF,2);

	if nargin < 4
		gpuFlag = false;
	end

	if nCol > 1
		x = zeros(NM, nCol);
		for i=1:nCol
			x(:,i) = ifftNM(XF(:,i),N,M);
		end
	else
		if gpuFlag
			XF = gpuArray(reshape(XF,N,M));
			x = transp(gather(ifft(XF(:, [K+1:M, 1:K]).')))*M;
		else
			XF = reshape(XF,N,M);
			x = transp(ifft(transp(XF(:, [K+1:M, 1:K]))))*M;
		end
		x = x(:);
	end
	
end
