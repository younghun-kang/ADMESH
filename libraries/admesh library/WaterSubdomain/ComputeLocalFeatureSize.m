function lfs = ComputeLocalFeatureSize(MA,D,dx)
%% ========================================================================
% Separate 1D & 2D area based on local feature size
%==========================================================================

if ~isempty(MA.BranchNodes)
id_MA = vertcat(MA.BranchNodes{:});
if islogical(id_MA)
    id_MA = find(id_MA);
end
else
    id_MA = [];
end

% bw_MA = zeros(MA.Size);
% bw_MA(id_MA) = 1;

[N,M] = size(D);
[I,J] = find(D);
[IMA,JMA] = ind2sub([N,M],id_MA);
K = sub2ind([N,M],I,J);

I = uint32(I);
J = uint32(J);
K = uint32(K);

%--------------------------------------------------------------------------
% Simpler and more memory-efficient... not exact
%--------------------------------------------------------------------------
MA2 = false(size(D));
MA2(id_MA) = 1;
d2MA = double(bwdist(MA2))*dx;

D1 = abs(D(K));
d2MA1 = abs(d2MA(K));
lfs1 = D1 + d2MA1;

lfs = sparse(I,J,lfs1,N,M);

%--------------------------------------------------------------------------
% Simpler way... not exact
%--------------------------------------------------------------------------
% [~,dist] = knnsearch([IMA,JMA],[I,J]);
% dist = dist*dx;
% lfs = sparse(I,J,abs(dist) + abs(D(K)),N,M); clear dist;

%--------------------------------------------------------------------------
% More accurate but too slow..
%--------------------------------------------------------------------------
% IN = full(D < 0);
% lfs = bwdistgeodesic(IN,id_MA);
% lfs(isinf(lfs)) = 0;
% 
% lfs = sparse(I,J,double(lfs(K)*dx) + abs(D(K)),N,M);
