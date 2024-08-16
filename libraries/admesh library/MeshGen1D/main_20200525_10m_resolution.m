%% Extract Open-channels using TopoToolbox
addpath(genpath('topotoolbox-master'));
path = 'U:\PR-NSF-PREEVENTS\DG-SAKE other watersheds\Bosque River Watershed\ArcGIS';
file = [path,'\ned10m_MBRW2.tif'];

DEM = GRIDobj(file);
DEMf = fillsinks(DEM);

FD = FLOWobj(DEMf);
A = flowacc(FD);

% figure; imageschs(DEM,min(gradient8(DEM),1));
% figure; imageschs(DEM,dilate(sqrt(A),ones(5)),'colormap',flipud(copper));

W = A > 1e5;
S = STREAMobj(FD,W);
S = klargestconncomps(S);

figure; imagesc(DEM); demcmap(DEM.Z); colorbar;
hold on; plot(S,'b','linewidth',1);

%% Read watershed boundary
path = 'U:\PR-NSF-PREEVENTS\DG-SAKE other watersheds\Bosque River Watershed\ArcGIS';
file = [path,'\MBRW_boundary_m.shp'];
MBRW = shaperead(file);

Bx = MBRW.X;
By = MBRW.Y;
BP = [Bx(:),By(:)];
BP(end,:) = []; % end point is NaN

%% Read 30m-DEM file in matrix format
% path = 'U:\PR-NSF-PREEVENTS\DG-SAKE other watersheds\Bosque River Watershed\ArcGIS';
% file = [path,'\ned30m31097_Clip.tif'];
% 
% [X,Y,Z] = ReadTifFile(file);
% 
% [sx,sy] = STREAMobj2XY(S);
% figure(1); 
% surf(X,Y,Z,'EdgeColor','none','FaceColor','interp');
% view(2); axis equal tight; mycolormap('dem'); colorbar;
% hold on; myPlotOnTop(sx,sy,'b');
% % myFancyPlot;
% 
% myPlotOnTop(Bx,By,'k','linewidth',1);

%% Read 10m-DEM file in matrix format
path = 'U:\PR-NSF-PREEVENTS\DG-SAKE other watersheds\Bosque River Watershed\ArcGIS';
file = [path,'\ned10m_MBRW2.tif'];
[X10,Y10,Z10] = ReadTifFile(file);
% remove out of area of interest
temp = X10(1,:) < 614457;
X10(:,temp) = []; Y10(:,temp) = []; Z10(:,temp) = [];
temp = X10(1,:) > 664588;
X10(:,temp) = []; Y10(:,temp) = []; Z10(:,temp) = [];
temp = Y10(:,1) < 3478036;
X10(temp,:) = []; Y10(temp,:) = []; Z10(temp,:) = [];
temp = Y10(:,1) > 3501227;
X10(temp,:) = []; Y10(temp,:) = []; Z10(temp,:) = [];

%--------------------------------------------------------------------------
% Plot
%--------------------------------------------------------------------------
figure;
surf(X10,Y10,Z10,'EdgeColor','none','FaceColor','interp');
view(2); axis equal tight; mycolormap('dem'); colorbar;


%% Re-arrange open-channels
[sx,sy] = STREAMobj2XY(S);
nanS = [0; find(isnan(sx))];
Check = cell(length(nanS)-1,1);
k = 0;
for i = 1 : length(nanS)-1
    x = sx(nanS(i)+1:nanS(i+1)-1);
    y = sy(nanS(i)+1:nanS(i+1)-1);

    [inS,onS] = inpolygon(x,y,Bx,By);
    if any(inS)
        Check{i} = find(inS);
        k = k + 1;     
        Sx{k} = x(inS);
        Sy{k} = y(inS);
        SP{k} = [Sx{k},Sy{k}];
    end
end

SX = []; SY = [];
for i = 1 : length(Sx)
    SX = vertcat(SX, Sx{i},nan);
    SY = vertcat(SY, Sy{i},nan);
end

figure; hold on;
for i = 1 : length(SP)
    myPlotOnTop(SP{i}(:,1),SP{i}(:,2),'b','linewidth',1);
end


%%
NewSP = CoarsenFlowline(SP,100);
NewBP = {CoarsenLineSegment(BP,1000)};

NewSX = []; NewSY = [];
for i = 1 : length(NewSP)
    NewSX = vertcat(NewSX, NewSP{i}(:,1),nan);
    NewSY = vertcat(NewSY, NewSP{i}(:,2),nan);
end


% figure; hold on;
% for i = 1 : length(NewSP)
%     myPlotOnTop(NewSP{i}(:,1),NewSP{i}(:,2),'b');
% end
% hold on; myPlotOnTop(NewBP{1}(:,1),NewBP{1}(:,2));


%% Read NHD_Flowline shapefile
path = 'U:\PR-NSF-PREEVENTS\DG-SAKE other watersheds\Bosque River Watershed\ArcGIS';
file = [path,'\NHD_Flowline_m.shp'];
NHDfl = shaperead(file);

figure; hold on;
surf(X,Y,Z,'EdgeColor','none','FaceColor','interp');
view(2); axis equal tight; mycolormap('dem'); colorbar;
myPlotOnTop(Bx,By,'k');
for i = 1 : length(NHDfl)
%     plot3(NHDfl(i).X,NHDfl(i).Y,1e3*ones(size(NHDfl(i).X)),'m')
    myPlotOnTop(NHDfl(i).X,NHDfl(i).Y,'m');
    
%     NHD_fl{i} = [NHDfl(i).X(:),NHDfl(i).Y(:)];
end












