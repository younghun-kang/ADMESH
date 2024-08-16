%==========================================================================
% main_MedialAxis_multiple_shorelines_v2
% - This function does this, this, and this...
% 
% Update history
% 2021-03-16 (v1) - written by Younghun Kang
% 2021-03-03 (v2) - clean up codes and change variable names
% 
%==========================================================================

%% ========================================================================
% Load datasets
%==========================================================================
clear all; clc; close all;
% addpath(genpath('admesh-lib-part'));
%--------------------------------------------------------------------------
% Load example file (variable "Points")
%--------------------------------------------------------------------------
ShorelineFile = 'D:\academic\data\gis\lower_neches\Shoreline_InBoundary.shp';
Shoreline     = shaperead(ShorelineFile);
% WBD = shaperead('D:\academic\data\gis\lower_neches\Boundary.shp');
% WBD_pgon = polyshape(WBD.X,WBD.Y);

for i = 1 : length(Shoreline)
    Shoreline(i).np = length(Shoreline(i).X);
end

%% ========================================================================
% Parameter setting
%==========================================================================
close all;

%%========================================================================
% Set boundary points from shorelines and set up parameters
%==========================================================================
fid = 0 : length(Shoreline)-1;
% fid = [1272,420,421,436,445,514,515,516,537,543,544,546,547,550,554,597];
% fid = [338, 340, 343, 345, 351, 352, 353, 358, 364, 366, 367, 381, 383, 386, 398, 400, 410, 412, 415, 417, 418, 419, 426, 440, 447, 449, 462, 470, 483, 484, 486, 490, 491, 496, 521, 524, 534, 540, 542, 551, 553, 556, 558, 559, 561, 1284];
% fid = 1571;
id = fid+1;

%--------------------------------------------------------------------------
% Construct polyshape with filtering small islands
%--------------------------------------------------------------------------
Shoreline1 = struct2table(Shoreline);
x = Shoreline1.X(id);
y = Shoreline1.Y(id);
x = cellfun(@(x) x(~isnan(x)),x,'UniformOutput',0);
y = cellfun(@(x) x(~isnan(x)),y,'UniformOutput',0);
A = cellfun(@(x,y) polyarea(x,y),x,y,'UniformOutput',0);
A = cell2mat(A);
I = A > 2000/(deg2km(1)*1e3)^2;
% I = A > 0;
pgon = polyshape(x(I),y(I));
pgon_removed = polyshape(x(~I),y(~I));

%--------------------------------------------------------------------------
% Load polyshape data (if desired)
%--------------------------------------------------------------------------
pgon_noholes = rmholes(pgon);

I = find(isnan(pgon.Vertices(:,1)));
I = [0; I(:); size(pgon.Vertices,1)+1];
XY_b = {pgon_noholes.Vertices}; % boundary including the dd-box
k1 = 1;
    
for i = 1 : length(I)-1
    i1 = I(i) + 1;
    i2 = I(i+1) - 1;
    J = [i1:i2, i1];
    
    if nnz(~ismember(XY_b{1},pgon.Vertices(J,:))) > 0
        k1 = k1 + 1;
        XY_b{k1} = pgon.Vertices(J,:);
    end
end
% x_b = XY_b(:,1); y_b = XY_b(:,2);

DEG2KM = 6378*pi/180;
m2deg = km2deg(1e-3);
% delta = 1e-5;
dx = 1*m2deg; % resolution of background grid (convert meter to degree)
hmin = 100*m2deg; % meter to degree
hmax = 10*hmin;

delta_cw = 100*m2deg; 
delta_lfs = delta_cw/2;
delta_Area = 1000*m2deg^2; % square meter to square degree (?) 
delta_A = delta_Area/dx^2; % conversion to number of pixels
pruning_rho = 4*delta_lfs;
pruning_theta = pi*.9;

%% ========================================================================
% Decompose domain
%==========================================================================
Nx = 10;
Ny = 10;
[dd_pgons,dd_Boxes,dd_ID,xg,yg] = DecomposePolyshape(pgon,Nx,Ny,dx,delta_cw*4);

% figure; hold on;
% for i = 1 : length(dd_pgons)
%     plot(dd_pgons(i).Vertices(:,1),dd_pgons(i).Vertices(:,2));
% %     plot(dd_Boxes(i));
% %     
% %     I = [dd_ID{i,2}; dd_ID{i,2}(end)*ones(size(dd_ID{i,1})); flip(dd_ID{i,2}); dd_ID{i,2}(1)*ones(size(dd_ID{i,1}))];
% %     J = [dd_ID{i,1}(1)*ones(size(dd_ID{i,2})); dd_ID{i,1}; dd_ID{i,1}(end)*ones(size(dd_ID{i,2})); flip(dd_ID{i,1})];
% %     plot(xg(I(:)),yg(J(:)),'r');
% end
% axis equal;


%% ========================================================================
% Compute VDT in DDM
%==========================================================================
[Vxg,Vyg,MaskW,MaskL] = ComputeVDT_DDM(pgon,dd_ID,dd_Boxes,xg,yg,dx);


%% ========================================================================
% Compute Medial Axis
%==========================================================================
MaskMA_raw = ComputeMA_DDM(Vxg,Vyg,dd_ID);

% [I,J] = find(MaskMA_raw);
% figure; plot(pgon.Vertices(:,1),pgon.Vertices(:,2),'k'); axis equal;
% hold on; plot(xg(J),yg(I),'.');

%% ========================================================================
% Thinning MA
%%=========================================================================
MaskMA_thinned = ThinningMA(MaskMA_raw);
fprintf('Thinning is done\n');

% [I,J] = find(MaskMA_thinned);
% figure; plot(pgon.Vertices(:,1),pgon.Vertices(:,2),'k'); axis equal;
% hold on; plot(xg(J),yg(I),'.');

%% ========================================================================
% Construct branches
%==========================================================================
MA_branch = ConstructMedialAxisBranch(MaskMA_thinned);
MA_branch = RemoveJointDuplicates(MA_branch);
MA = MA_branch;
fprintf('Branch construction is done\n');

%% ========================================================================
% Pruning corners
%==========================================================================
MA_pruned = SkeletonPruning_Choi_level1branch_v2(Vxg,Vyg,MA,pruning_rho,pruning_theta,dx);
fprintf('Pruning corner branches is done\n');

%% ========================================================================
% Compute distance map to medial axis
%==========================================================================
MA = MA_pruned;
D2MA = DistanceToMA(MA,{MaskL,MaskW},dd_ID,dx);

%% ========================================================================
% Compute lfs
%==========================================================================
lfs = sparse(size(MaskW,1),size(MaskW,2));
wbar = waitbar(0);
for iDD = 1 : size(dd_ID,1)
    I = dd_ID{iDD,1};
    J = dd_ID{iDD,2};
  
    lfs(I,J) = sqrt(Vxg(I,J).^2 + Vyg(I,J).^2) + abs(D2MA(I,J));
    
    waitbar(iDD/size(dd_ID,1),wbar,sprintf('Wairbar (%d/%d)',iDD,size(dd_ID,1)));
end
delete(wbar);
fprintf('Computing local feature size is done\n');
clear Dg Vxg Vyg;
% clear D2MA;
%% ========================================================================
% Split 1D and 2D masks based on lfs
%==========================================================================
[Water2D,Water1D] = SplitMask(MaskW,lfs,delta_lfs,delta_A,MA_pruned);
[Land2D,Land1D] = SplitMask(MaskL,lfs,delta_lfs,delta_A,MA_pruned);

% figure; imshow(Water2D,'XData',xg,'YData',yg); set(gca,'YDir','Normal');
% figure; imshow(Water1D,'XData',xg,'YData',yg); set(gca,'YDir','Normal');
% figure; imshow(Land2D,'XData',xg,'YData',yg); set(gca,'YDir','Normal');
% figure; imshow(Land1D,'XData',xg,'YData',yg); set(gca,'YDir','Normal');

%% ========================================================================
% Fill 2D mask regions
%==========================================================================
delta_filling = delta_lfs*sqrt(2);
[Mask2DCell,Mask1DCell] = FillMask2D({Water2D,Land2D},{Water1D,Land1D},MA_pruned,lfs,delta_filling,delta_A,xg,yg,XY_b,dx);
% [Mask2DCell,Mask1DCell] = FillMask2D(Land2D,Land1D,MA_pruned,lfs,delta_filling,delta_A,xg,yg,XY_b,delta);

% figure; imshow(flipud(Mask2DCell{1}));
% figure; imshow(flipud(Mask1DCell{1}));
% figure; imshow(flipud(Mask2DCell{2}));
% figure; imshow(flipud(Mask1DCell{2}));

Water2D = Mask2DCell{1};
Water1D = Mask1DCell{1};
Land2D = Mask2DCell{2};
Land1D = Mask1DCell{2};

clear Mask2DCell Mask1DCell;


%% ========================================================================
% Transfer 1D masks to 2D masks based on connectivity.
%==========================================================================
% CCW1D = bwconncomp(Water1D);
% CCL1D = bwconncomp(Land1D);

Ratio = .7;

BranchNodesID = vertcat(MA_pruned.BranchNodes{:});
BranchNodesID = unique(BranchNodesID);

CC = bwconncomp(Land1D);
if CC.NumObjects > 0
%--------------------------------------------------------------------------
% Compute PA ratio
%--------------------------------------------------------------------------
CCRP = regionprops(CC,'MaxFeretProperties','MinFeretProperties',...
    'Perimeter','Area','ConvexArea','BoundingBox','SubarrayIdx',...
    'MajorAxisLength','MinorAxisLength','Centroid');
CCRP = struct2table(CCRP,'AsArray',1);
CC.PAR = CCRP.Perimeter./CCRP.Area; % Isoperimetric ratio
CC.IPR = (CCRP.Perimeter).^2./CCRP.Area; % Isoperimetric ratio
CC.FDR = CCRP.MaxFeretDiameter./CCRP.MinFeretDiameter; % Feret diameter ratio
CC.CAR = CCRP.ConvexArea./CCRP.Area;
CC.AR = CCRP.MajorAxisLength./CCRP.MinorAxisLength;


% %--------------------------------------------------------------------------
% % Use square buffer with size 100m and check ratio of water2d and land2d
% %--------------------------------------------------------------------------
% BoundingBox = CCRP.BoundingBox;
% BoundingBox = mat2cell(BoundingBox,ones(size(BoundingBox,1),1),size(BoundingBox,2));
% 
% temp_I = cellfun(@(x) floor(x(2))-100:floor(x(2))+x(4)+100,BoundingBox,'UniformOutput',0);
% temp_J = cellfun(@(x) floor(x(1))-100:floor(x(1))+x(3)+100,BoundingBox,'UniformOutput',0);
% 
% temp_I = cellfun(@(x) x(x >= 1 & x <= CC.ImageSize(1)),temp_I,'UniformOutput',0);
% temp_J = cellfun(@(x) x(x >= 1 & x <= CC.ImageSize(2)),temp_J,'UniformOutput',0);
% 
% temp_WLratio = cellfun(@(x,y) nnz(Water2D(x,y))/nnz(Land2D(x,y)),temp_I,temp_J);
% temp_WLratio(isnan(temp_WLratio)) = 10;
% 
% temp_K = temp_WLratio > 1 & CC.IPR > 30;

%--------------------------------------------------------------------------
% Use Euclidean buffer with size MinorAxisLength and check ratio of water2d and land2d
%--------------------------------------------------------------------------
Buffer = cell(CC.NumObjects,1);
for i = 1 : CC.NumObjects
    id = CC.PixelIdxList{i};
    [temp_I,temp_J] = ind2sub(MA_pruned.Size,id);
    
    Box = CCRP.BoundingBox(i,:);
%     BufferSize = round(CCRP.MaxFeretDiameter(i)/2);
%     BufferSize = round(max(Box(3:4))/2);
    BufferSize = round(CCRP.MinorAxisLength(i)/2);
    
    minI = max(1,floor(Box(2))-BufferSize);
    maxI = min(MA_pruned.Size(1),floor(Box(2)+Box(4))+BufferSize);
    Buffer_I = minI : maxI;
    
    minJ = max(1,floor(Box(1))-BufferSize);
    maxJ = min(MA_pruned.Size(2),floor(Box(1)+Box(3))+BufferSize);
    Buffer_J = minJ : maxJ;
    
    [Buffer_I,Buffer_J] = meshgrid(Buffer_I,Buffer_J);
    Buffer_I = Buffer_I(:);
    Buffer_J = Buffer_J(:);
    
    [~,dist] = knnsearch([temp_I,temp_J],[Buffer_I(:),Buffer_J(:)]);
    
    I = dist < BufferSize & dist > 0;
    if nnz(I) == 0
        0;
    end
    Buffer{i} = sub2ind(MA_pruned.Size,Buffer_I(I),Buffer_J(I));
end
temp_WLratio = cellfun(@(x) nnz(Water2D(x))/nnz(Land2D(x)),Buffer);
I = cellfun(@(x) isempty(x),Buffer);
temp_WLratio(I) = 10;
temp_WLratio(isnan(temp_WLratio)) = 0;

temp_K = temp_WLratio > 2 & CC.IPR > 30;

% maxD2MA = cellfun(@(x) full(max(D2MA(x))),CC.PixelIdxList);
% meanD2MA = cellfun(@(x) full(mean(D2MA(x))),CC.PixelIdxList);
% temp_K = temp_WLratio > 2 & maxD2MA(:) < delta_cw*0.5;% & CC.IPR > 30;

temp_K = vertcat(CC.PixelIdxList{temp_K});
L1DtoW2D = false(CC.ImageSize(1),CC.ImageSize(2));
L1DtoW2D(temp_K) = 1;
% figure; imshow(L1DtoW2D,'XData',xg,'YData',yg); set(gca,'YDir','Normal');

else
    L1DtoW2D = false(CC.ImageSize(1),CC.ImageSize(2));
end
    
%%



i = 3;    
figure; hold on;
id = Buffer{i};
[temp_I,temp_J] = ind2sub(MA_pruned.Size,id);
plot(xg(temp_J),yg(temp_I),'.');

id = CC.PixelIdxList{i};
[temp_I,temp_J] = ind2sub(MA_pruned.Size,id);
plot(xg(temp_J),yg(temp_I),'.');


temp_bufferXY = [];
for i = 1 : length(Buffer)
temp = false(MA_pruned.Size);
temp(Buffer{i}) = 1;
temp = Mask2XY(temp,xg,yg);
temp_bufferXY = [temp_bufferXY; temp{1}; nan nan];
end

figure; plot(temp_bufferXY(:,1),temp_bufferXY(:,2),'k');
%% ========================================================================
% Construct new water mask
%%=========================================================================
%--------------------------------------------------------------------------
% Filter out isolated regions
%--------------------------------------------------------------------------
NewW2D = Water2D | L1DtoW2D;
CC = bwconncomp(NewW2D);
for i = 1 : CC.NumObjects
    id = CC.PixelIdxList{i};
    if ~any(Water2D(id))
        L1DtoW2D(id) = 0;
    end
end
clear NewW2D;

%--------------------------------------------------------------------------
% Identify branch nodes indexes in L1DtoW2D mask
%--------------------------------------------------------------------------
BranchNodesID = vertcat(MA_pruned.BranchNodes{:});
BranchNodesID = unique(BranchNodesID);
K = BranchNodesID(L1DtoW2D(BranchNodesID));

L1DtoW2Dskel = false(size(L1DtoW2D));
L1DtoW2Dskel(K) = 1;

MA_LandBarrier = ConstructMedialAxisBranch(L1DtoW2Dskel);
MA_LandBarrier = RemoveJointDuplicates(MA_LandBarrier);
% MA_LandBarrier = PruneLevel1BranchByLength(dx,MA_LandBarrier,10/(DEG2KM*1e3));
% MA_LandBarrier = MergeShortBridge(MA_LandBarrier,10/(DEG2KM*1e3),dx);
% MA_LandBarrier = RemoveJointDuplicates(MA_LandBarrier);

% clear L1DtoW2Dskel;

%--------------------------------------------------------------------------
% Construct new water mask
%--------------------------------------------------------------------------
% NewW2D = Water2D | Water1D | L1DtoW2D;
% NewW2D = Water2D | L1DtoW2D;

%% ========================================================================
% Patch Water 2D mask with L1DtoW2D mask
%==========================================================================
L1DtoW2D_XY = Mask2XY(L1DtoW2D,xg,yg);
if ~isempty(L1DtoW2D_XY)
pgon_x = cellfun(@(x) x(:,1),L1DtoW2D_XY,'UniformOutput',0);
pgon_y = cellfun(@(x) x(:,2),L1DtoW2D_XY,'UniformOutput',0);
pgon_x = cellfun(@(x) [x; nan],pgon_x,'UniformOutput',0);
pgon_y = cellfun(@(x) [x; nan],pgon_y,'UniformOutput',0);

MiterLimit = 2; % This value is tentative. Experiments would be need
pgon_add = polyshape(pgon_x,pgon_y);
pgon_add = polybuffer(pgon_add,dx*sqrt(2),'JointType','miter','MiterLimit',MiterLimit);

pgon_new = union(pgon,pgon_add);
else
    pgon_new = pgon;
end
% Fill small holes
pgon_new_holes = holes(pgon_new);
pgon_new_holes_area = area(pgon_new_holes);
I = pgon_new_holes_area/(dx^2) < 10;
if nnz(I) > 0
pgon_new_holes = union(pgon_new_holes(I));
pgon_new = union(pgon_new,pgon_new_holes);
end
% Update polygon (backup old one)
pgon_old = pgon;
pgon = pgon_new;

%% ========================================================================
% Begin again with updated polygon
%==========================================================================
pgon_noholes = rmholes(pgon);

I = find(isnan(pgon.Vertices(:,1)));
I = [0; I(:); size(pgon.Vertices,1)+1];
XY_b = {pgon_noholes.Vertices}; % boundary including the dd-box
k1 = 1;
    
for i = 1 : length(I)-1
    i1 = I(i) + 1;
    i2 = I(i+1) - 1;
    J = [i1:i2, i1];
    
    if nnz(~ismember(XY_b{1},pgon.Vertices(J,:))) > 0
        k1 = k1 + 1;
        XY_b{k1} = pgon.Vertices(J,:);
    end
end


%% Decompose domain
[dd_pgons,dd_Boxes,dd_ID,xg,yg] = DecomposePolyshape(pgon,Nx,Ny,dx,delta_cw*4);

%% ========================================================================
% Compute VDT in DDM
%==========================================================================
[Vxg,Vyg,MaskW,MaskL] = ComputeVDT_DDM(pgon,dd_ID,dd_Boxes,xg,yg,dx);
% Land mask is not required at this step
Vxg(MaskL) = 0;
Vyg(MaskL) = 0;

%% Compute Medial Axis
MaskMA_raw = ComputeMA_DDM(Vxg,Vyg,dd_ID);

% [I,J] = find(MaskMA_raw);
% figure; plot(pgon.Vertices(:,1),pgon.Vertices(:,2),'k'); axis equal;
% hold on; plot(xg(J),yg(I),'.');

%% ========================================================================
% Thinning MA
%%=========================================================================
MaskMA_thinned = ThinningMA(MaskMA_raw);
fprintf('Thinning is done\n');

%% ========================================================================
% Construct branches (with redundant joints)
%==========================================================================
MA_branch = ConstructMedialAxisBranch(MaskMA_thinned);
MA_branch = RemoveJointDuplicates(MA_branch);
MA = MA_branch;
fprintf('Branch construction is done\n');

%% ========================================================================
% Pruning corners
%==========================================================================
MA_pruned = SkeletonPruning_Choi_level1branch_v2(Vxg,Vyg,MA,pruning_rho,pruning_theta,dx);
fprintf('Pruning corner branches is done\n');

%% ========================================================================
% Compute distance map to medial axis
%==========================================================================
MA = MA_pruned;
D2MA = DistanceToMA(MA,{MaskW},dd_ID,dx);

%% ========================================================================
% Compute lfs
%==========================================================================
lfs = sparse(size(MaskW,1),size(MaskW,2));
wbar = waitbar(0);
for iDD = 1 : size(dd_ID,1)
    I = dd_ID{iDD,1};
    J = dd_ID{iDD,2};
  
    lfs(I,J) = sqrt(Vxg(I,J).^2 + Vyg(I,J).^2) + abs(D2MA(I,J));

    waitbar(iDD/size(dd_ID,1),wbar,sprintf('Wairbar (%d/%d)',iDD,size(dd_ID,1)));
end
delete(wbar);
fprintf('Computing local feature size is done\n');
clear Dg D2MA Vxg Vyg;

%% ========================================================================
% Find water masks based on lfs (filter with connectivity to MA)
% (This is a part of script below, but written in function)
%==========================================================================
[Water2D,Water1D] = SplitMask(MaskW,lfs,delta_lfs,delta_A,MA_pruned);

%% ========================================================================
% Fill 2D mask regions
%==========================================================================
tempW2D = Water2D;
tempW2D(L1DtoW2D) = 1;
tempW1D = Water1D;
tempW1D(L1DtoW2D) = 1;

temp_MA = MA_pruned;
temp_MA.BranchNodes = vertcat(temp_MA.BranchNodes(:),MA_LandBarrier.BranchNodes(:));

delta_filling = delta_lfs*sqrt(2);
[~,~,M1DtoM2DCell] = FillMask2D({tempW2D},{tempW1D},temp_MA,lfs,delta_filling,delta_A,xg,yg,XY_b,dx);

M1DtoM2D = M1DtoM2DCell{1};

tempW2D = Water2D;
tempW2D(M1DtoM2D) = 1;
tempW1D = Water1D;
tempW1D(M1DtoM2D) = 0;

BranchNodesID = vertcat(temp_MA.BranchNodes{:});
BranchNodesID = unique(BranchNodesID);

[I,J] = ind2sub(size(tempW2D),BranchNodesID);
K = tempW1D(BranchNodesID);
MA1D = sparse(I,J,K,size(tempW2D,1),size(tempW2D,2));

[I,J] = ind2sub(size(tempW2D),BranchNodesID);
K = tempW2D(BranchNodesID);
MA2D = sparse(I,J,K,size(tempW2D,1),size(tempW2D,2));

CC = bwconncomp(full(tempW1D));
wbar = waitbar(0);
fwbar = @(x,y) waitbar(x/y,wbar,sprintf('Wairbar (%d/%d)',x,y));
for i = 1 : CC.NumObjects
    id = CC.PixelIdxList{i};
    if ~any(MA1D(id))
        tempW2D(id) = 1;
        tempW1D(id) = 0;
    end
    fwbar(i,CC.NumObjects);
end
delete(wbar);

%----------------------------------------------------------------------
% Remove 2D mask not containing any MA points (it may need only when 
% the filling results isolated mask, which should not happen)
%----------------------------------------------------------------------
CC = bwconncomp(full(tempW2D));
for i = 1 : CC.NumObjects
    id = CC.PixelIdxList{i};
    if ~any(MA2D(id)) || numel(id) < delta_A
        tempW2D(id) = 0;
        tempW1D(id) = 1;
    end
end

Water2D = tempW2D;
Water1D = tempW1D;

%% ========================================================================
% Construct MA1D
%==========================================================================
MA_1D = ConstructMA1D(MA_pruned,Water1D,Water2D);
MA_1D = RemoveJointDuplicates(MA_1D);
% MA_1D = PruneLevel1BranchByLength(dx,MA_1D,10*m2deg);
% MA_1D = MergeShortBridge(MA_1D,10*m2deg,dx);
% MA_1D = RemoveJointDuplicates(MA_1D);

% ! WIP: Pruning L1 branch with < 10 length isn't good.. it loose some branch connected to W2D
% ! The reason added this was that very short branch connec
%% ========================================================================
% Extract Water2D mask boundary
%==========================================================================
[~,W2D_boundaryID] = Mask2XY(Water2D,xg,yg);

temp_x = cellfun(@(x) xg(x(:,2)),W2D_boundaryID,'UniformOutput',0);
temp_y = cellfun(@(x) yg(x(:,1)),W2D_boundaryID,'UniformOutput',0);
pgon_W2D = polyshape(temp_x,temp_y);
%% ========================================================================
% Connect branches to Water2D boundary
%==========================================================================
MA = MA_1D;
MA_connected = ConnectMA1Dto2DArea_v2(MA,W2D_boundaryID,2);
MA_connected = PruneLevel1BranchByLength(dx,MA_connected,100/(deg2km(1)*1e3),MA_connected.FlagConnected2D);
MA_connected = RemoveJointDuplicates(MA_connected);
% MA_connected = ConnectMA1Dto2DArea_v2(MA_connected,W2D_boundaryID,hmin/2/dx);

MA = MA_LandBarrier;
if ~isempty(MA_LandBarrier.BranchNodes)
MA_connectedLand = ConnectMA1Dto2DArea_v2(MA,W2D_boundaryID,2);
MA_connectedLand = PruneLevel1BranchByLength(dx,MA_connectedLand,100/(deg2km(1)*1e3),MA_connectedLand.FlagConnected2D);
MA_connectedLand = RemoveJointDuplicates(MA_connectedLand);
else
    MA_connectedLand = MA_LandBarrier;
end
% MA_connectedLand = ConnectMA1Dto2DArea_v2(MA_connectedLand,W2D_boundaryID,hmin/2/dx);

%% ========================================================================
% Construct mainstreams
%==========================================================================
addpath('D:\academic\PR-PREEVENTS\codes_MATLAB\ADMESH\ADMESH_v4\libraries\MeshGen1D');

MA = MA_connected;
MA_1D_mainstream = ConstructMainStream_v3(MA);
MA_1D_mainstream = RemoveJointDuplicates(MA_1D_mainstream);

MA = MA_connectedLand;
MA_1D_mainstreamLand = ConstructMainStream_v3(MA);
MA_1D_mainstreamLand = RemoveJointDuplicates(MA_1D_mainstreamLand);
rmpath('D:\academic\PR-PREEVENTS\codes_MATLAB\ADMESH\ADMESH_v4\libraries\MeshGen1D');


%% ========================================================================
% Construct projection maps for 1D domain and 2D boundary
%==========================================================================
addpath('D:\academic\PR-PREEVENTS\codes_MATLAB\ADMESH\MeshGen1D');
Settings1D.h_min          = hmin;
Settings1D.h_max          = hmax;
Settings1D.h0             = 1*hmin;
Settings1D.K              = 20;
Settings1D.g              = 0.15;
Settings1D.SmoothingRMSE = 5*km2deg(1e-3); % set negative value to disable smoothing

I = [vertcat(MA_1D_mainstream.BranchNodes{:});
    vertcat(MA_1D_mainstreamLand.BranchNodes{:})];
% I =  vertcat(MA_LandBarrier.BranchNodes{:});
[~,J] = unique(I);
J = setdiff(1:length(I),J);
I = I(J);
[I,J] = ind2sub(MA_1D_mainstream.Size,I);
fixedPoints1 = [xg(J),yg(I)];

%----------------------------------------------------------------------
% Add fixed points for junctions
%----------------------------------------------------------------------
f = @(x) ind2xy(x.BranchNodes,x.Size,xg,yg);
MA_1D_mainstream.BranchXY = f(MA_1D_mainstream);
MA_1D_mainstreamLand.BranchXY = f(MA_1D_mainstreamLand);
W2D_boundary = sub2xy(W2D_boundaryID,xg,yg);

POINTS = {MA_1D_mainstream.BranchXY;
    MA_1D_mainstreamLand.BranchXY;
    W2D_boundary};

% POINTS1 = vertcat(POINTS{:});
% POINTS1 = vertcat(POINTS1{:});
% POINTS2 = unique(POINTS1,'rows');
% fixedPoints = [];
% nfixedP = 0;
% for i = 1 : size(POINTS2,1)
%     if nnz(ismember(POINTS1,POINTS2(i,:),'rows')) > 1
%         nfixedP = nfixedP + 1;
%         fixedPoints(nfixedP,:) = POINTS2(i,:);
%     end
% end

% id = find(MA.BranchLength < Setting1D.h_min);
% Points(id) = [];

fixedPoints2 = [];
nfixedP = 0;
for i = 1 : 2
    Points = POINTS{i};
    for j = 1 : length(Points)
        nfixedP = nfixedP + 1;
        fixedPoints2(nfixedP,:) = [Points{j}(1,1),Points{j}(1,2)];
        nfixedP = nfixedP + 1;
        fixedPoints2(nfixedP,:) = [Points{j}(end,1),Points{j}(end,2)];
    end
end
fixedPoints = unique([fixedPoints1; fixedPoints2],'rows');

clear PI0;
for j = 1 : length(POINTS)
k = 0;
wbar = waitbar(0);
fwbar = @(x,y) waitbar(x/y,wbar,sprintf('Wairbar (%d/%d)',x,y));
Points = POINTS{j};
for i = 1 : length(Points)
    x = Points{i}(:,1);
    y = Points{i}(:,2);
    
    PI1 = ComputePathCurvature([x,y],fixedPoints,Settings1D.SmoothingRMSE,'Warning','Off');
%     PI(i) = PI1;
    if ~isempty(PI1)
        k = k + 1;
        PI0{j}(k) = PI1;
    end
    fwbar(i,length(Points));
end
delete(wbar);
end

figure; hold on; daspect([1 1 1]);
C = {'b',mycolors('bn'),mycolors('dgn')};
for j = 1 : length(PI0)
    PI1 = PI0{j};
    for i = 1 : length(PI1)
        plot(PI1(i).x1(PI1(i).p),PI1(i).y1(PI1(i).p),'-o','color',C{j});
%         myPlot3(PI1(i).x1(PI1(i).p),PI1(i).y1(PI1(i).p),PI1(i).k(PI1(i).p));
%         myPlot3(PI1(i).sx(PI1(i).p),PI1(i).sy(PI1(i).p),PI1(i).k(PI1(i).p));
    end
end

% figure; hold on; daspect([1 1 1]);
% PI1 = PI0{3};
% for i = 1 : length(PI1)
% %     plot(PI1(i).x1(PI1(i).p),PI1(i).y1(PI1(i).p));
% %     myPlot3(PI1(i).x1(PI1(i).p),PI1(i).y1(PI1(i).p),PI1(i).k(PI1(i).p));
%     myPlot3(PI1(i).sx(PI1(i).p),PI1(i).sy(PI1(i).p),PI1(i).k(PI1(i).p));
% end


%% ========================================================================
% Construct dummy constraints for ADMESH
%==========================================================================
ADMESH_folder = 'D:\academic\PR-PREEVENTS\codes_MATLAB\ADMESH\ADMESH_v4\';
Filename = [ADMESH_folder,'test_files\test\simple_example'];
%----------------------------------------------------------------------
% Setup projection to use
%----------------------------------------------------------------------
for j = 1 : length(PI0)
    for i = 1 : length(PI0{j})
        PI0{j}(i).x = PI0{j}(i).x1;
        PI0{j}(i).y = PI0{j}(i).y1;
    end
end
% PI1 = PI(1:length(PI1));
% PI2 = PI(length(PI1)+1:end);

%----------------------------------------------------------------------
% Setup ADMESH input
%----------------------------------------------------------------------
PTS = [];

xl = min(pgon.Vertices(:,1)) - Settings1D.h_min;
xr = max(pgon.Vertices(:,1)) + Settings1D.h_min;
yl = min(pgon.Vertices(:,2)) - Settings1D.h_min;
yr = max(pgon.Vertices(:,2)) + Settings1D.h_min;
external_boundary = {[xl, yl; xr, yl; xr, yr; xl, yr; xl, yl]};
for i = 1 : length(external_boundary)
    PTS.Poly(i).x = external_boundary{i}(:,1);
    PTS.Poly(i).y = external_boundary{i}(:,2);
end

k = 0;
for j = 1 : length(PI0)
    PI = PI0{j};
    for i = 1 : length(PI)
        p = PI(i).p;
        x1 = PI(i).x(p);
        y1 = PI(i).y(p);
        k1 = PI(i).k(p);
        
        k = k + 1;
        PTS.Constraints(k).num = 18;
        PTS.Constraints(k).xy = [x1,y1];
        PTS.Constraints(k).type = 'line';
        PTS.Constraints(k).data = [];
        PTS.Constraints(k).Kappa = k1;
    end
end

Settings = [];
Settings.Res = 'Low';
Settings.hmax = Settings1D.h_max;
Settings.hmin = Settings1D.h_min;
Settings.K.Status = 'on';
Settings.K.Value = Settings1D.K;

Settings.G.Status = 'on';
Settings.G.Value = Settings1D.g/sqrt(1);
Settings.View.Status = 'on';

Settings.R.Status = 'off';
Settings.B.Status = 'off';
Settings.T.Status = 'off';
Settings.DummyConstraint = 1;

DummyConstFile = [Filename,'_dummyconstraint'];
folder = fileparts(DummyConstFile);
if ~exist(folder,'dir')
    mkdir(folder);
end
save(DummyConstFile,'PTS','Settings');

cd(ADMESH_folder);

%% ========================================================================
% Mesh generation 1D with updated mesh size
%==========================================================================
% Points = MA.XY(52);

OverwriteFile = [Filename,'_dummyconstraint_updated_1d_hsize'];

load(OverwriteFile,'X','Y','h');
figure; hold on; daspect([1 1 1]);
for i = 1 : length(PI0)
    PI = PI0{i};
    if isempty(PI)
        continue;
    end
temp = struct2table(PI,'AsArray',1);
temp = table2struct(temp,'ToScalar',1);
temp_x = cellfun(@(x,y) [x(y);nan],temp.x,temp.p,'UniformOutput',0);
temp_y = cellfun(@(x,y) [x(y);nan],temp.y,temp.p,'UniformOutput',0);
temp_x = vertcat(temp_x{:});
temp_y = vertcat(temp_y{:});
temp_h = interp2(X,Y,h,temp_x,temp_y);
myPlot3(temp_x,temp_y,temp_h);
end


%----------------------------------------------------------------------
% Add fixed points for junctions
%----------------------------------------------------------------------
fixedPoints = [];
nfixedP = 0;
for i = 1 : length(PI0)
    PI = PI0{i};
    for j = 1 : length(PI)
        nfixedP = nfixedP + 1;
        x = PI(j).x1(PI(j).p);
        y = PI(j).y1(PI(j).p);
        k = PI(j).k(PI(j).p);
        fixedPoints(nfixedP,:) = [x(1), y(1)];
        fixedPointsK(nfixedP) = k(1);
        nfixedP = nfixedP + 1;
        fixedPoints(nfixedP,:) = [x(end), y(end)];
        fixedPointsK(nfixedP) = k(end);
    end
end

clear Mesh1D0;

wbar = waitbar(0);
fwbar = @(x,y) waitbar(x/y,wbar,sprintf('Wairbar (%d/%d)',x,y));
for j = 1 : length(PI0)
    PI1 = PI0{j};
    for i = 1 : length(PI1)
        x = PI1(i).x1(PI1(i).p);
        y = PI1(i).y1(PI1(i).p);
        
        id = ismember([x,y],fixedPoints,'rows');
        fixedPoints1D = PI1(i).p(id);
        
        %----------------------------------------------------------------------
        % Generate 1D mesh
        %----------------------------------------------------------------------
        mesh1d = MeshGeneration1D(PI1(i),fixedPoints1D,Settings1D,OverwriteFile);
        Mesh1D0{j}(i) = mesh1d;
        fwbar(i,length(PI1));
    end
end
delete(wbar);


figure; hold on; axis equal;
C = {'b','r','m'};
for i = 1 : length(Mesh1D0)
    PI1 = PI0{i};
    if isempty(PI1)
        continue;
    end
    temp = struct2table(PI1,'AsArray',1);
    temp = table2struct(temp,'ToScalar',1);
    temp_x = cellfun(@(x,y) [x(y);nan],temp.x,temp.p,'UniformOutput',0);
    temp_y = cellfun(@(x,y) [x(y);nan],temp.y,temp.p,'UniformOutput',0);
    temp_x = vertcat(temp_x{:});
    temp_y = vertcat(temp_y{:});
    plot(temp_x,temp_y,'k');
    
    Mesh1D1 = Mesh1D0{i};
    temp = struct2table(Mesh1D1);
    temp = table2struct(temp,'ToScalar',1);
    temp_x = cellfun(@(x) [x;nan],temp.X,'UniformOutput',0);
    temp_y = cellfun(@(x) [x;nan],temp.Y,'UniformOutput',0);
    temp_x = vertcat(temp_x{:});
    temp_y = vertcat(temp_y{:});
    plot(temp_x,temp_y,'-o','color',C{i},'MarkerFaceColor',C{i},'MarkerSize',4);
end
% 
figure; hold on; axis equal;
Mesh1D1 = Mesh1D0{1};
for i = 1 : length(Mesh1D1)
    plot(Mesh1D1(i).X,Mesh1D1(i).Y,'-*');
end
% axis([-94.0173  -94.0006   30.0604   30.0735]);

%% ========================================================================
% Postprocess 0: Transfer shoreline constraints to internal boundaries if it is straightline
%==========================================================================
temp_Mesh1D = Mesh1D0;
temp_Mesh1D = cellfun(@(x) x(:),temp_Mesh1D,'UniformOutput',0);
Mesh1Dpost0 = [];

Mesh1D1 = Mesh1D0{3};
I = [];
for i = 1 : length(Mesh1D1)
    if isequal([Mesh1D1(i).X(1),Mesh1D1(i).Y(1)],[Mesh1D1(i).X(end),Mesh1D1(i).Y(end)])...
        && length(unique([Mesh1D1(i).X,Mesh1D1(i).Y],'rows')) == 2
        Mesh1D1(i).X = [];
        Mesh1D1(i).Y = [];
        I = [I;i];
    end
end

Mesh1Dpost0 = temp_Mesh1D;
Mesh1Dpost0{2} = [Mesh1Dpost0{2}; temp_Mesh1D{3}(I)];
Mesh1Dpost0{3}(I) = [];
%% ========================================================================
% Postprocess 1: Merge 1D point clusters which have neighboring fixed 
% points within hmin/2
%==========================================================================
temp_Mesh1D = Mesh1Dpost0;
Mesh1Dpost1 = [];

Points1D = [];
for i = 1 : length(temp_Mesh1D)
    Mesh1D1 = temp_Mesh1D{i};
    if isempty(Mesh1D1)
        continue;
    end
    Mesh1D1 = struct2table(Mesh1D1,'AsArray',1);
    Points1D = [Points1D; vertcat(Mesh1D1.X{:}), vertcat(Mesh1D1.Y{:})];
    
%     temp_x = cellfun(@(x) x([1 end]),Mesh1D1.X,'UniformOutput',0);
%     temp_y = cellfun(@(x) x([1 end]),Mesh1D1.Y,'UniformOutput',0);
%     Points1D = [Points1D; vertcat(temp_x{:}), vertcat(temp_y{:})];
end
Points1D = unique(Points1D,'rows');

[id,dist] = knnsearch(Points1D,Points1D,'k',5);
id = id(:,2:end);
dist = dist(:,2:end);
id1 = [];
for i = 1 : size(dist,1)
    k = find(dist(i,:) < hmin/4 & dist(i,:) > 0,1,'last');
    id1{i} = sort([i; id(i,1:k)']);
end
id1 = id1(:);
I = cellfun(@(x) length(x) > 1,id1);
id1 = id1(I);

for i = 1 : length(id1)
    for j = i+1 : length(id1)
        if all(ismember(id1{j},id1{i}))
            id1{j} = [];
        end
        if all(ismember(id1{i},id1{j}))
            id1{i} = [];
        end
    end
end
I = cellfun(@(x) ~isempty(x),id1); 
id1 = id1(I);

for i = 1 : length(id1)
    for j = i+1 : length(id1)
        if any(ismember(id1{j},id1{i}))
            error('Something is wrong.');
        end
    end
end

mfp = [];
iMap = [];
for i = 1 : length(id1)
    mfp(i,:) = mean(Points1D(id1{i},:));
    iMap = [iMap; i*ones(length(id1{i}),1)];
end
mfp(:,3) = 0;

iFixedMerge = vertcat(id1{:});
for j = 1 : length(temp_Mesh1D)
    Mesh1D1 = temp_Mesh1D{j};
    if isempty(Mesh1D1)
        continue;
    end
    for i = 1 : length(Mesh1D1)
        x = Mesh1D1(i).X;
        y = Mesh1D1(i).Y;
        K = Mesh1D1(i).K;
        
        [I,J] = ismember([x,y],Points1D(iFixedMerge,:),'rows');
        
%         k = false(length(x),1);
%         k = [];
        if any(I)
            I = find(I);
            J = nonzeros(J);
            J = iMap(J);
            
            x(I) = mfp(J,1);
            y(I) = mfp(J,2);
            K(I) = mfp(J,3);
%             k(I) = 1;
            k1 = find(diff(I) == 1);
            k2 = I(k1);
            k = k2(J(k1) == J(k1+1));
            x(k) = [];
            y(k) = [];
            K(k) = [];
        end
        Mesh1D1(i).X = x;
        Mesh1D1(i).Y = y;
        Mesh1D1(i).K = K;
    end
    
    Mesh1D1 = struct2table(Mesh1D1,'AsArray',1);
    I = cellfun(@(x) length(x),Mesh1D1.X) > 1;
    Mesh1D1 = Mesh1D1(I,:);
    Mesh1D1 = table2struct(Mesh1D1);
    
    Mesh1Dpost1{j} = Mesh1D1;
end

%% ========================================================================
% Post-process 2: Remove 1D elements nodes with length smaller than hmin/2
%==========================================================================
temp_Mesh1D = Mesh1Dpost1;
Mesh1Dpost2 = [];

fixedPoints = [];
nfixedP = 0;
for i = 1 : length(temp_Mesh1D)
    Mesh1D1 = temp_Mesh1D{i};
    for j = 1 : length(Mesh1D1)
        nfixedP = nfixedP + 1;
        x = Mesh1D1(j).X;
        y = Mesh1D1(j).Y;
        k = Mesh1D1(j).K;
        fixedPoints(nfixedP,:) = [x(1), y(1)];
        fixedPointsK(nfixedP) = k(1);
        nfixedP = nfixedP + 1;
        fixedPoints(nfixedP,:) = [x(end), y(end)];
        fixedPointsK(nfixedP) = k(end);
    end
end

for j = 1 : length(temp_Mesh1D)
    Mesh1D1 = temp_Mesh1D{j};
    for i = 1 : length(Mesh1D1)
        x = Mesh1D1(i).X;
        y = Mesh1D1(i).Y;
        k = Mesh1D1(i).K;
        while 1
            d = sqrt(diff(x).^2 + diff(y).^2);
            I = find(d < hmin/2);
            %         [~,I] = min(d);
            %         if d(I) > hmin*10
            %             break;
            %         end
            iFixed = find(ismember([x,y],fixedPoints,'rows'));
            I = setdiff(I,iFixed);
            if isempty(I)
                break;
            end
            [~,J] = min(d(I));
            I = I(J);
            x(I(1)) = [];
            y(I(1)) = [];
            k(I(1)) = [];
        end
        Mesh1Dpost2{j}(i).X = x;
        Mesh1Dpost2{j}(i).Y = y;
        Mesh1Dpost2{j}(i).K = k;
    end
    
end

% %% ========================================================================
% % Postprocess 2: Merge 1d elements which have length smaller than hmin/2 
% % and having fixed points on both end
% %==========================================================================
% temp_Mesh1D = Mesh1D0;
% Mesh1Dpost2 = [];
% 
% [id,dist] = knnsearch(fixedPoints,fixedPoints,'k',2);
% id = id(:,end);
% dist = dist(:,end);
% I = dist < hmin/4 & dist > 0;
% I1 = find(I);
% I2 = id(I1);
% 
% mfp = (fixedPoints(I1,:) + fixedPoints(I2,:))/2;
% mfpk = (fixedPointsK(I1) + fixedPointsK(I2))/2;
% % hold on; plot(fixedPoints(I,1),fixedPoints(I,2),'oc');
% % hold on; plot(mfp(:,1),mfp(:,2),'*c');
% 
% 
% iFixedMerge = [I1(:); I2(:)];
% for j = 1 : length(temp_Mesh1D)
%     Mesh1D1 = temp_Mesh1D{j};
%     for i = 1 : length(Mesh1D1)
%         x = Mesh1D1(i).X;
%         y = Mesh1D1(i).Y;
%         K = Mesh1D1(i).K;
%         
%         [I,J] = ismember([x,y],fixedPoints(I1,:),'rows');
%         
% %         k = false(length(x),1);
% %         k = [];
%         if any(I)
%             I = find(I);
%             J = nonzeros(J);
%             
%             x(I) = mfp(J,1);
%             y(I) = mfp(J,2);
%             K(I) = mfpk(J);
% %             k(I) = 1;
%             k = I(diff(I) == 1);
%             x(k) = [];
%             y(k) = [];
%             K(k) = [];
%         end
%         
%         [I,J] = ismember([x,y],fixedPoints(I2,:),'rows');
%         
%         if any(I)
%             I = find(I);
%             J = nonzeros(J);
%             
%             x(I) = mfp(J,1);
%             y(I) = mfp(J,2);
%             K(I) = mfpk(J);
% %             k(I) = 1;
%             
%             k = I(diff(I) == 1);
%             x(k) = [];
%             y(k) = [];
%             K(k) = [];
%         end
%         
% %         x(k) = [];
% %         y(k) = [];
% %         K(k) = [];
%         
%         Mesh1D1(i).X = x;
%         Mesh1D1(i).Y = y;
%         Mesh1D1(i).K = K;
%     end
%     
%     Mesh1D1 = struct2table(Mesh1D1);
%     I = cellfun(@(x) length(x),Mesh1D1.X) > 2;
%     Mesh1D1 = Mesh1D1(I,:);
%     Mesh1D1 = table2struct(Mesh1D1);
%     
%     Mesh1Dpost2{j} = Mesh1D1;
% end
% 
%     
% figure; hold on; axis equal;
% C = {'-*b','-*r','-*m'};
% for i = 1 : length(Mesh1Dpost2)
%     PI1 = PI0{i};
%     temp = struct2table(PI1);
%     temp = table2struct(temp,'ToScalar',1);
%     temp_x = cellfun(@(x,y) [x(y);nan],temp.x,temp.p,'UniformOutput',0);
%     temp_y = cellfun(@(x,y) [x(y);nan],temp.y,temp.p,'UniformOutput',0);
%     temp_x = vertcat(temp_x{:});
%     temp_y = vertcat(temp_y{:});
%     plot(temp_x,temp_y,'k');
%     
%     Mesh1D1 = Mesh1Dpost2{i};
%     temp = struct2table(Mesh1D1);
%     temp = table2struct(temp,'ToScalar',1);
%     temp_x = cellfun(@(x) [x;nan],temp.X,'UniformOutput',0);
%     temp_y = cellfun(@(x) [x;nan],temp.Y,'UniformOutput',0);
%     temp_k = cellfun(@(x) [x;nan],temp.K,'UniformOutput',0);
%     temp_x = vertcat(temp_x{:});
%     temp_y = vertcat(temp_y{:});
%     temp_k = vertcat(temp_k{:});
%     plot(temp_x,temp_y,C{i});
% %     myScatter3(temp_x,temp_y,temp_k);
% end
% hold on; plot(fixedPoints([I1;I2],1),fixedPoints([I1;I2],2),'oc');
% hold on; plot(mfp(:,1),mfp(:,2),'*c');


%% ========================================================================
% Postprocess 3: Remove 1D constraints that completely included in shoreline constraints
%==========================================================================
temp_Mesh1D = Mesh1Dpost2;
Mesh1Dpost3 = [];

temp = temp_Mesh1D{3};
temp = struct2table(temp);
if iscell(temp)
temp = [vertcat(temp.X{:}),vertcat(temp.Y{:})];
else
    temp = [temp.X,temp.Y];
end

for j = 1 : 2
    Mesh1D1 = temp_Mesh1D{j};
    if isempty(Mesh1D1)
        continue;
    end
    for i = 1 : length(Mesh1D1)
        x = Mesh1D1(i).X;
        y = Mesh1D1(i).Y;
        K = Mesh1D1(i).K;
        
        [~,dist] = knnsearch(temp,[x,y]);
%         if all(ismember([x,y],temp,'rows'))
%         if nnz(~ismember([x,y],temp,'rows')) < 3
        
%         I = dist < hmin/2;
%         I([1 end]) = 0;
%         Mesh1D1(i).X(I) = [];
%         Mesh1D1(i).Y(I) = [];
%         Mesh1D1(i).K(I) = [];
        if max(dist) < hmin/2
            Mesh1D1(i).X = [];
            Mesh1D1(i).Y = [];
            Mesh1D1(i).K = [];
        end
    end
    Mesh1D1 = struct2table(Mesh1D1);
    if iscell(Mesh1D1.X)
    I = cellfun(@(x) length(x),Mesh1D1.X) > 1;
    Mesh1D1 = Mesh1D1(I,:);
    Mesh1D1 = table2struct(Mesh1D1);
    else
        Mesh1D1 = table2struct(Mesh1D1,'ToScalar',1);
    end
    
    Mesh1Dpost3{j} = Mesh1D1;
end
Mesh1Dpost3{3} = temp_Mesh1D{3};

%%
figure; hold on; axis equal;
C = {'b',mycolors('bn'),mycolors('dgn')};
% plot(pgon_W2D);
for i = 1 : length(Mesh1Dpost3)
    PI1 = PI0{i};
    if isempty(PI1)
        continue;
    end
    temp = struct2table(PI1,'AsArray',1);
    temp = table2struct(temp,'ToScalar',1);
    if iscell(temp.x)
    temp_x = cellfun(@(x,y) [x(y);nan],temp.x,temp.p,'UniformOutput',0);
    temp_y = cellfun(@(x,y) [x(y);nan],temp.y,temp.p,'UniformOutput',0);
    temp_x = vertcat(temp_x{:});
    temp_y = vertcat(temp_y{:});
    else
        temp_x = temp.x;
        temp_y = temp.y;
    end
%     plot(temp_x,temp_y,'k');
    
    Mesh1D1 = Mesh1Dpost3{i};
    temp = struct2table(Mesh1D1);
    temp = table2struct(temp,'ToScalar',1);
    if iscell(temp.X)
    temp_x = cellfun(@(x) [x;nan],temp.X,'UniformOutput',0);
    temp_y = cellfun(@(x) [x;nan],temp.Y,'UniformOutput',0);
    temp_k = cellfun(@(x) [x;nan],temp.K,'UniformOutput',0);
    temp_x = vertcat(temp_x{:});
    temp_y = vertcat(temp_y{:});
    temp_k = vertcat(temp_k{:});
    else
        temp_x = temp.X;
        temp_y = temp.Y;
        temp_k = temp.K;
    end
    plot(temp_x,temp_y,'-o','color',C{i},'MarkerFaceColor',C{i},'MarkerSize',4);
%     myScatter3(temp_x,temp_y,temp_k);
end

% figure; hold on; axis equal;
% Mesh1D1 = Mesh1Dpost3{1};
% for i = 1 : length(Mesh1D1)
%     plot(Mesh1D1(i).X,Mesh1D1(i).Y,'-*');
% end
% axis([-94.0173  -94.0006   30.0604   30.0735]);

%% Write ADMESH input for single-line constraint (pass both 1D lines and 2D boundaries as constraints)
temp_Mesh_write = Mesh1Dpost3;

PTS = [];

xl = min(pgon.Vertices(:,1)) - Settings1D.h_min;
xr = max(pgon.Vertices(:,1)) + Settings1D.h_min;
yl = min(pgon.Vertices(:,2)) - Settings1D.h_min;
yr = max(pgon.Vertices(:,2)) + Settings1D.h_min;
external_boundary = {[xl, yl; xr, yl; xr, yr; xl, yr; xl, yl]};

for i = 1 : length(external_boundary)
    PTS.Poly(i).x = external_boundary{i}(:,1);
    PTS.Poly(i).y = external_boundary{i}(:,2);
end

k = 0;
ConstraintNum = [18 17 19];
for j = 1 : length(temp_Mesh_write)
    Mesh1D1 = temp_Mesh_write{j};
    for i = 1 : length(Mesh1D1)
        k = k + 1;
        PTS.Constraints(k).num = ConstraintNum(j);
        PTS.Constraints(k).xy = [Mesh1D1(i).X,Mesh1D1(i).Y];
        PTS.Constraints(k).type = 'line';
        PTS.Constraints(k).data = [];
        PTS.Constraints(k).Kappa = Mesh1D1(i).K;
    end
end


Settings = [];
Settings.Res = 'Low';
Settings.hmax = Settings1D.h_max;
Settings.hmin = Settings1D.h_min;
Settings.K.Status = 'on';
Settings.K.Value = Settings1D.K;

Settings.G.Status = 'on';
Settings.G.Value = Settings1D.g/sqrt(1);
Settings.View.Status = 'on';

Settings.R.Status = 'off';
Settings.B.Status = 'off';
Settings.T.Status = 'off';

Settings.DummyConstraint = 0;

cd(ADMESH_folder);
ConstraintFilename = [Filename,'_single_constraint'];
save(ConstraintFilename,'PTS','Settings','Settings1D');











