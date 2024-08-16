function AreaInfo = SeparateAreas(MA,D,dx,delta_lfs,delta_A)
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

bw_MA = zeros(MA.Size);
bw_MA(id_MA) = 1;

dMA = bwdist(bw_MA)*dx;

id_out = D > 0 | isnan(D);
dMA(id_out) = nan;

dMA = double(dMA);
lfs = (abs(dMA) + abs(D));
%--------------------------------------------------------------------------
% Truncate based on local feature size
%--------------------------------------------------------------------------
lfs_2D = lfs;
lfs_2D(id_out) = nan;
id_L1D = lfs_2D < delta_lfs;
lfs_2D(id_L1D) = nan;

%--------------------------------------------------------------------------
% Get boundary
%--------------------------------------------------------------------------
bw_2D = zeros(MA.Size);
bw_2D(~id_L1D) = 1;
bw_2D(id_out) = 0;
[B,L2D] = bwboundaries(bw_2D); % without 'noholes' option, bwboundaries result islands as new layers, which means islands become water bodies in our case
id_island = bw_2D == 0 & L2D > 0;
[~,id_island] = bwboundaries(id_island); % fill holes in the id_island
id_island = id_island > 0;
L2D(id_island) = 0; % no needed if 'noholes' option is used above
L1D(id_island) = 0;

% bw_2D = lfs_2D;
% bw_2D(~isnan(bw_2D)) = 1;
% bw_2D(isnan(bw_2D)) = 0;
% bw_2D = bw_2D';
% [B,L2D] = bwboundaries(bw_2D);
% L2D(bw_2D == 0) = 0;
% L2D = L2D';

% L2D(bw_2D == 0) = 0;
%--------------------------------------------------------------------------
% Compute area of each labeled region
%--------------------------------------------------------------------------
nL = zeros(max(L2D(:)),1);
for i = 1 : max(L2D(:))
    nL(i) = nnz(L2D == i);
end
A = nL;

%--------------------------------------------------------------------------
% Remove small regions
%--------------------------------------------------------------------------
id = find(A < delta_A);
B(id) = [];
for i = 1 : length(id)
    L2D(L2D == id(i)) = 0;
end
%--------------------------------------------------------------------------
% Re-label 2D layers
%--------------------------------------------------------------------------
% [L2D_boundary,L2D] = bwboundaries(L2D);
idL2D = unique(L2D);
k = 0;
for i = 1 : length(idL2D)
    id = find(L2D == idL2D(i));
    if ~isempty(id)
        L2D(id) = k;
        k = k + 1;
    end
end

%--------------------------------------------------------------------------
% Make 1D layer
%--------------------------------------------------------------------------
L0 = double(~id_out);
L1D = L0;
L1D = full(L1D);
L1D(L2D > 0) = 0;
L1D(id_island) = 0;

%--------------------------------------------------------------------------
% Re-label 1D layers
%--------------------------------------------------------------------------
% [L1D_boundary,L1D] = bwboundaries(L1D);
idL1D = unique(L1D);
k = 0;
for i = 1 : length(idL1D)
    id = find(L1D == idL1D(i));
    if ~isempty(id)
        L1D(id) = k;
        k = k + 1;
    end
end


%--------------------------------------------------------------------------
% Get indices for MA1D and MA2D
%--------------------------------------------------------------------------
id_MA1D = intersect(id_MA,find(L1D > 0));
id_MA2D = intersect(id_MA,find(L2D > 0));
% XY_MA1D = [X(id_MA1D),Y(id_MA1D)];
% XY_MA2D = [X(id_MA2D),Y(id_MA2D)];

AreaInfo.id_MA1D = id_MA1D;
AreaInfo.L1D     = L1D;
% AreaInfo.XY_MA1D = XY_MA1D;

AreaInfo.id_MA2D = id_MA2D;
AreaInfo.L2D     = L2D;
% AreaInfo.XY_MA2D = XY_MA2D;

AreaInfo.Island = id_island;



