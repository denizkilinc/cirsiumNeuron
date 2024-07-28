function XF = fftNM(x, N, M, gpuFlag)
    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------
    
	% N... number of variables
	% M... number of time points
	K = (M-1)/2; % number of freqs
	NM = N*M;
	nCol = size(x,2);

	if nargin < 4
		gpuFlag = false;
	end

	if nCol > 1
		XF = zeros(NM, nCol);
		for i=1:nCol
			XF(:,i) = fftNM(x(:,i),N,M);
		end
	else
		if gpuFlag == true
			x = gpuArray(reshape(x,N,M));
			XF = transp(gather(fft(x.')))/M;
		else
			x = full(reshape(x,N,M));
			XF = transp(fft(transp(x)))/M;
		end
		XF = XF(:, [K+2:M, 1:K+1]);
		XF = XF(:);
	end
	
end
