%% ========================================================================
% Use watershed obtained from existing mesh
%==========================================================================

PTS = [];
SHP = shaperead('test_files\GCEW\boundary.shp');
PTS.Poly.x = SHP.X;
PTS.Poly.y = SHP.Y;

Fort14 = Read14File('test_files\GCEW\Goodwin_Creek_Watershed_Mesh_90_g13.14');
[~,P] = freeBoundary(Fort14.TR);
PTS.Poly.x = P(:,1);
PTS.Poly.y = P(:,2);

Fort14 = Read14File('test_files\GCEW\Goodwin_Creek_Watershed_Mesh_90_g13.14');
load('test_files\GCEW\Channel_Width_Interpolant.mat');
F = scatteredInterpolant(Fw.Points(:,1),Fw.Points(:,2),Fw.Values,'nearest','nearest');

P = Fort14.TR.Points;
Constraints = Fort14.Constraints;
for i = 1 : length(Constraints)
    x = P(Constraints(i).nodeStr,1);
    y = P(Constraints(i).nodeStr,2);
    Constraints(i).data(:,1) = 1;
    Constraints(i).data(:,2) = F(x,y);
    Constraints(i).data(:,3) = F(x,y);

    Constraints(i).num = -18;
    Constraints(i).xy = [x(:), y(:)];
end

PTS.Constraints = Constraints;

save('test_files\GCEW\GCEW2','PTS');

%% ========================================================================
% Use watershed boundary extracted from DEM with TopoToolbox
%==========================================================================
dem_file = 'D:\academic\data\gis\goodwin_creek\dem\USGS_1_n35w090_NAD83_UTM16N_clip3.tif';
%--------------------------------------------------------------------------
% Find watershed boundary with TopoToolbox
%--------------------------------------------------------------------------
DEM = GRIDobj(dem_file);
FD  = FLOWobj(DEM);
S   = STREAMobj(FD,'minarea',1000);
S   = klargestconncomps(S,1);
D = drainagebasins(FD,S);
[sx,sy] = STREAMobj2XY(S);
[x,y] = getcoordinates(D);
[X,Y] = meshgrid(x,y);
I = (D.Z == 1);

%--------------------------------------------------------------------------
% Modify watershed boundary with buffer along open channels
%--------------------------------------------------------------------------
[~,dist] = knnsearch([sx,sy],[X(:),Y(:)]);
I_ex = dist < 200;
I_ex = reshape(I_ex,size(X,1),size(X,2));
I = I | I_ex;

%--------------------------------------------------------------------------
% Get new watershed boundary
%--------------------------------------------------------------------------
x1 = X(I);
y1 = Y(I);
k = boundary(x1,y1,1);
x1 = x1(k);
y1 = y1(k);

%--------------------------------------------------------------------------
% Convert open channel variables to struct
%--------------------------------------------------------------------------
FL = NaNdlm2struct([sx,sy],'Boundary',[x1(:),y1(:)]);

%--------------------------------------------------------------------------
% Compute watershed centroid to shift open-channel data and gauge locations
%--------------------------------------------------------------------------
pgon1 = polyshape(x1,y1);
Fort14 = Read14File('test_files\GCEW\Goodwin_Creek_Watershed_Mesh_90_g13.14');
[~,P] = freeBoundary(Fort14.TR);
pgon2 = polyshape(P(:,1),P(:,2));
[xc1,yc1] = centroid(pgon1);
[xc2,yc2] = centroid(pgon2);
dx = xc1 - xc2;
dy = yc1 - yc2;
load('test_files\GCEW\Channel_Width_Interpolant.mat');
F = scatteredInterpolant(Fw.Points(:,1),Fw.Points(:,2),Fw.Values,'nearest','nearest');
F1 = scatteredInterpolant(Fw.Points(:,1)+dx,Fw.Points(:,2)+dy,Fw.Values,'nearest','nearest');

%--------------------------------------------------------------------------
% Set new ADMESH input file
%--------------------------------------------------------------------------
PTS = [];
PTS.Poly.x = x1;
PTS.Poly.y = y1;
Fort14 = Read14File('test_files\GCEW\Goodwin_Creek_Watershed_Mesh_90_g13.14');
load('test_files\GCEW\Channel_Width_Interpolant.mat');
F = scatteredInterpolant(Fw.Points(:,1),Fw.Points(:,2),Fw.Values,'nearest','nearest');
Constraints = [];
Constraints(length(FL),1).num = [];
for i = 1 : length(Constraints)
    Constraints(i).num = -18;
    Constraints(i).xy = FL{i};
    Constraints(i).data(:,1) = ones(size(FL{i},1),1);
    Constraints(i).data(:,2) = F1(FL{i});
    Constraints(i).data(:,3) = F1(FL{i});
end
PTS.Constraints = Constraints;
save('test_files\GCEW\GCEW4','PTS');

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
figure; hold on; axis equal;
imageschs(DEM,[],'colormap',[1 1 1],'colorbar',false)
h = plotdbfringe(FD,S,'colormap',parula,'width',30);
plot(S,'b')
% figure;
plot(x1,y1,'k','linewidth',2); axis equal;
% PTS = app.PTS;
% hold on; plot(PTS.Poly.x,PTS.Poly.y,'r','linewidth',2);
% figure; hold on; axis equal;
plot(pgon1);
plot(pgon2);
plot(xc1,yc1,'*b');
plot(xc2,yc2,'*r');
dx = xc1 - xc2;
dy = yc1 - yc2;
plot(RO_Stations.xs+dx,RO_Stations.ys+dy,'ok','markerFaceColor','w');
% for i = 1 : length(Constraints)
% x = P(Constraints(i).nodeStr,1);
% y = P(Constraints(i).nodeStr,2);
% plot(x,y,'r');
% end




