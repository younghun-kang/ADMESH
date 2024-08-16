function [MA_1D,AreaInfo_Filled] = Fill2DArea_v2(MA,AreaInfo,Vmag,dx)
%% ========================================================================
% Fill 2D area
%==========================================================================
id_MA1D = AreaInfo.id_MA1D;
L1D     = AreaInfo.L1D;
                          
id_MA2D = AreaInfo.id_MA2D;
L2D     = AreaInfo.L2D;

id_island = AreaInfo.Island;

id_out = full(Vmag == 0);

L2D_Filled = L2D;
L1D_Filled = L1D;
for i = 1 : length(id_MA2D)
    k = id_MA2D(i);
    iL2D = L2D(k);
    
%     x = X(k); y = Y(k);
    iVmag = Vmag(k);
    
    num_nghb_buffer = ceil(iVmag/dx)+1;
    
    [I,J] = ind2sub(size(Vmag),k);
    
    I1 = max(I - num_nghb_buffer,1);
    I2 = min(I + num_nghb_buffer,size(Vmag,1));
    
    J1 = max(J - num_nghb_buffer,1);
    J2 = min(J + num_nghb_buffer,size(Vmag,2));
    
    I3 = I1:I2;
    J3 = J1:J2;
    [I3,J3] = meshgrid(I3,J3);
    I3 = I3(:);
    J3 = J3(:);
    
    k_nghb = sub2ind(size(Vmag),I3,J3);
    
%     k_nghb = k + [-N-1, -N, -N+1, -1, 1, N-1, N, N+1];
    
%     tempx = X(k_nghb);
%     tempy = Y(k_nghb);
    
%     dist = sqrt( (tempx - x).^2 + (tempy - y).^2);
    dist = dx*sqrt((J3 - J).^2 + (I3 - I).^2);
    id = (dist < iVmag);
    id = k_nghb(id);
    id(id_out(id)) = [];
%     I4 = I3(id);
%     J4 = J3(id);
    
    L2D_Filled(id) = iL2D;
    L1D_Filled(id) = 0;
end
%--------------------------------------------------------------------------
% Re-label layers as some layers can be merged
%--------------------------------------------------------------------------
k1 = 0;
k2 = 0;
for i = 1 : max(max(L1D_Filled(:)),max(L2D_Filled(:)))
    id = L1D_Filled == i;
    if nnz(id) > 0
        k1 = k1 + 1;
        L1D_Filled(id) = k1;
    end
    
    id = L2D_Filled == i;
    if nnz(id) > 0
        k2 = k2 + 1;
        L2D_Filled(id) = k2;
    end
end
L2D_Filled(L2D_Filled > 0) = 1;
L1D_Filled(L1D_Filled > 0) = 1;

[~,L2D_Filled] = bwboundaries(L2D_Filled);
[~,L1D_Filled] = bwboundaries(L1D_Filled);

L2D_Filled(id_island) = 0;
L1D_Filled(id_island) = 0;
L2D_Filled(id_out) = 0;
L1D_Filled(id_out) = 0;

id_MA1D_Filled = intersect(id_MA1D, find(L1D_Filled > 0));
id_MA2D_Filled = intersect(id_MA2D, find(L2D_Filled > 0));

AreaInfo_Filled.L2D     = L2D_Filled;
AreaInfo_Filled.L1D     = L1D_Filled;

AreaInfo_Filled.id_MA1D = id_MA1D_Filled;
AreaInfo_Filled.id_MA2D = id_MA2D_Filled;

%% ========================================================================
% Remove MA points in 2D area
%==========================================================================
nBranch = 0;
BranchNodes = [];
JointConnectivity = [];
% XY = [];
for i = 1 : length(MA.BranchNodes)
    id = MA.BranchNodes{i};
    is_MA1D = L1D_Filled(id) > 0;
    if nnz(is_MA1D) > 1 % exclude single point
        nBranch = nBranch + 1;
        id = id(is_MA1D);
        BranchNodes{nBranch} = id;
        JointConnectivity(nBranch,:) = [id(1), id(end)];
%         XY{nBranch} = [X(id),Y(id)];
    end
end

MA_1D.nBranch = nBranch;
MA_1D.Size = MA.Size;
MA_1D.BranchNodes = BranchNodes;
MA_1D.JointConnectivity = JointConnectivity;
% MA_1D.XY = XY;



