addpath(genpath('D:\academic\PR-PREEVENTS\codes_MATLAB\MeshGeneration1D_v2'));

DEMFile = 'D:\academic\data\DEM\east_fork_white_river\3dep_13arc\Clip.tif';
BoundaryFile = 'D:\academic\data\gis\east_fork_white_river\AOI.shp';

FL = ExtractOpenChannels(DEMFile,1e4);

DEM = ReadTifFile(DEMFile);
fDEM = griddedInterpolant(DEM.X',DEM.Y',DEM.Z');

[sx,sy] = STREAMobj2XY(FL);

figure; hold on;
surf(DEM.X,DEM.Y,DEM.Z,'EdgeColor','None'); view(2); axis equal;
myPlot3(sx,sy,'b');



sx = sx(~isnan(sx));
sy = sy(~isnan(sy));
[PARAMS,CS] = RetrieveCrossSectionParameters_v2(fDEM,[sx,sy],0.2e-2,1e-5);


%%
figure;
AXmesh = axes;
hold on;
x = [PARAMS.Txl; flipud(PARAMS.Bxl)];
y = [PARAMS.Tyl; flipud(PARAMS.Byl)];
z = [PARAMS.BEL + PARAMS.H; flipud(PARAMS.BEL)];
patch(x,y,z,'edgecolor','interp');

x = [flipud(PARAMS.Bxl); PARAMS.Bxr];
y = [flipud(PARAMS.Byl); PARAMS.Byr];
z = [flipud(PARAMS.BEL); PARAMS.BEL];
patch(x,y,z,'edgecolor','interp');

x = [PARAMS.Bxr; flipud(PARAMS.Txr)];
y = [PARAMS.Byr; flipud(PARAMS.Tyr)];
z = [PARAMS.BEL; flipud(PARAMS.BEL+PARAMS.H)];
patch(x,y,z,'edgecolor','interp');
mycmap('viridis'); colorbar;
axis equal;


Dist = ComputeFlowlineDistance(sx,sy);

figure; hold on;
plot(Dist,fDEM(sx,sy));
plot(Dist,PARAMS.BEL);
plot(Dist,PARAMS.Tyl+PARAMS.BEL)
plot(Dist,PARAMS.Tyr+PARAMS.BEL)