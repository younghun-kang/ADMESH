clear all; clc; close all;
%--------------------------------------------------------------------------
% Input arguments
%--------------------------------------------------------------------------
DEMFile = 'examples\MBRW_ned30m.tif';
BoundaryFile = 'examples\MBRW_boundary.shp';

Setting.h_min     = 30;
Setting.h_max     = 2e2;
Setting.h0        = 1*Setting.h_min;
Setting.K         = 20;
Setting.g         = 0.2;
Setting.Smoothing = false;

%% ========================================================================
% Read data and smoothing
%==========================================================================
%--------------------------------------------------------------------------
% Extract Open-channels using TopoToolbox
%--------------------------------------------------------------------------
FL = ExtractOpenChannels(DEMFile,1e4);

%--------------------------------------------------------------------------
% Read watershed boundary
%--------------------------------------------------------------------------
BP = shaperead(BoundaryFile);
BP = [BP.X(:), BP.Y(:)];

%--------------------------------------------------------------------------
% Smoothing using TopoToolbox (CRS algorithm)
%--------------------------------------------------------------------------
k = 5;
SmoothFL = SmoothingCRS(FL,BP,k);

%% ========================================================================
% Generate mesh of flowlines with CRS smoothing
%==========================================================================
Mesh1D = [];
for i = 1 : length(SmoothFL)
    x = SmoothFL{i}(:,1);
    y = SmoothFL{i}(:,2);
    
    %----------------------------------------------------------------------
    % Add fixed points for junctions
    %----------------------------------------------------------------------
    fixedPoints = [];
    for j = 1 : length(SmoothFL)
        fixedPoints(j,:) = [SmoothFL{j}(end,1),SmoothFL{j}(end,2)];
    end
    id = ismember(fixedPoints,[x,y],'rows');
    fixedPoints = fixedPoints(id,:);
    
    %----------------------------------------------------------------------
    % Generate 1D mesh
    %----------------------------------------------------------------------
    mesh1d = MeshGeneration1D([x,y],fixedPoints,Setting);
    Mesh1D.X{i,1} = mesh1d.X(:);
    Mesh1D.Y{i,1} = mesh1d.Y(:);
    Mesh1D.K{i,1} = mesh1d.Kappa(:);

%     mesh1d = MeshGeneration1D_v2([x,y],fixedPoints,Setting,'Middle_Bosque_river_watershed_updated_1d_hsize');
%     Mesh1D.X{i,1} = mesh1d.X(:);
%     Mesh1D.Y{i,1} = mesh1d.Y(:);
%     Mesh1D.K{i,1} = mesh1d.K(:);
%     Mesh1D.PI{i,1} = mesh1d.PI;

end
%% ========================================================================
% Plot results
%==========================================================================
figure; hold on;
plot(BP(:,1),BP(:,2),'k','linewidth',2);
% plot(FL,'b','linewidth',2);
for i = 1 : length(SmoothFL)
%     plot(SmoothFL{i}(:,1),SmoothFL{i}(:,2),'k');
%     plot(Mesh1D.X{i},Mesh1D.Y{i},'-ob','MarkerFaceColor','b','linewidth',2);
    myPlot3(Mesh1D.X{i},Mesh1D.Y{i},1./(Setting.K*Mesh1D.K{i}),'linewidth',8);
end
daspect([1 1 1]);


%%

Channels = cell(length(Mesh1D.X),1);
s = linspace(0,1,1e3)';
for i = 1 : length(Channels)
    Channels{i} = [Mesh1D.X{i},Mesh1D.Y{i},Mesh1D.K{i}];
%     Channels{i} = [Mesh1D.PI{i}.x(s),Mesh1D.PI{i}.y(s),Mesh1D.PI{i}.k(s)];
end

% Channels = [];

Setting2D.Res = 'Low';
Setting2D.hmax = Setting.h_max;
Setting2D.hmin = Setting.h_min;
Setting2D.K.Status = 'on';
Setting2D.K.Value = Setting.K;

Setting2D.G.Status = 'on';
Setting2D.G.Value = Setting.g/sqrt(1);
Setting2D.View.Status = 'on';

Setting2D.R.Status = 'off';
Setting2D.B.Status = 'off';
Setting2D.T.Status = 'off';

cd('D:\academic\PR-PREEVENTS\codes_MATLAB\ADMESH_yk\ADMESH version 2.0_channel');
% cd('D:\academic\PR-PREEVENTS\codes_MATLAB\ADMESH_yk\ADMESH version 2.0_magnetic_force');
% WriteADMESH_input('test',BoundaryXY,Channels,Setting2D);
% WriteADMESH_input('test_step1_1dmesh',BoundaryXY,Channels,Setting2D);
WriteADMESH_input('Middle_Bosque_river_watershed',{BP},Channels,Setting2D);


