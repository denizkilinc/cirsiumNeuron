function [DN, DNinv] = dftmtxNM(N,M)
	% Generate extended DFT matrix.
	% This indexing scheme seems complicated but it is just handiwork. 
    % This is ~2 orders of magnitude faster than doing this with a loop.
    
    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------
    
	K = (M-1)/2;

	aIdx = (1:M*N)';
	aIdx = aIdx(:,ones(1,M));
	bIdx = reshape(1:M*N,N,M);
	bIdx = bIdx((1:N)'*ones(1,M),:);
	cIdx = (1:M^2);
	cIdx = cIdx(ones(1,N),:);
	%d = 1/M*exp(-1i*2*pi*(-K:K)'*fc*tVec);
	d = 1/M*exp(-1i*2*pi*(-K:K)'*(0:M-1)/M);

	DN = sparse(aIdx(:),bIdx(:),d(cIdx(:)));
	DNinv = M*DN';
end
