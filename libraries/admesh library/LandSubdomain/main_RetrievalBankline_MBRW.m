DEMFile = 'D:\academic\data\gis\middle_bosque_river\dem\usgs_3dep_13arc_Clip.tif';
% DEMFile = 'D:\academic\data\gis\middle_bosque_river\dem\3dep_13arc_Clip_ProjectRast.tif';

%% Read 10m-DEM file in matrix format
TIF = ReadTifFile(DEMFile);
X = TIF.X;
Y = TIF.Y;
Z = TIF.Z;

%%
% load Mesh1D_MBRW_GCS_NAD83;
load('D:\academic\PR-PREEVENTS\codes_MATLAB\ADMESH\ADMESH_v3\test_files\MBRW\GCS_NAD83\MBRW_GCA_NAD83_single_constraint.mat')
% load('D:\academic\PR-PREEVENTS\codes_MATLAB\ADMESH\ADMESH_v3\test_files\MBRW\GCS_NAD83\MBRW_GCA_NAD83_hmin200hmax1000_single_constraint.mat')

Mesh1D = [];
for i = 1 : length(PTS.Constraints)
    Mesh1D(i).X = PTS.Constraints(i).xy(:,1);
    Mesh1D(i).Y = PTS.Constraints(i).xy(:,2);
end
%%
m2deg = km2deg(1e-3);
fDEM = griddedInterpolant(X',Y',Z');

for i = 1 : length(Mesh1D)
cs_delta = 0.1*m2deg;
ChannelElem2DWidth = 30*m2deg;
Centerline = [Mesh1D(i).X,Mesh1D(i).Y];
% Centerline = [PI(i).x(PI(i).p),PI(i).y(PI(i).p)];
[PARAMS_bank{i}, CS_bankline{i}] = RetrieveBankline(fDEM,Centerline,ChannelElem2DWidth,cs_delta,2,'degree');
PARAMS_banksym{i} = RetrieveBankline(fDEM,Centerline,ChannelElem2DWidth,cs_delta,5,'degree');
% [PARAMS_banksym{i},CS_bankline{i}] = RetrieveBankline(fDEM,Centerline,ChannelElem2DWidth,cs_delta,6,'degree');
end

%%
i = 2;
figure;
hold on;
plot(PARAMS_bank{i}.x,PARAMS_bank{i}.y,'k');
temp_x = [PARAMS_bank{i}.Bxr; nan; PARAMS_bank{i}.Bxl];
temp_y = [PARAMS_bank{i}.Byr; nan; PARAMS_bank{i}.Byl];
plot(temp_x,temp_y,'b');
temp_x = [PARAMS_banksym{i}.Bxr; nan; PARAMS_banksym{i}.Bxl];
temp_y = [PARAMS_banksym{i}.Byr; nan; PARAMS_banksym{i}.Byl];
plot(temp_x,temp_y,'r');
axis equal;

%% Debug
deg2m = deg2km(1)*1e3;

XYq = cursor_info.Position([1 2]);
j = knnsearch([PARAMS_bank{i}.x,PARAMS_bank{i}.y],XYq);
% j = 1794;
% j = 431;


figure; hold on;
% plot(csd,csz,'k');
plot(CS_bankline{i}(j).d*deg2m,CS_bankline{i}(j).z,'k');
plot(PARAMS_bank{i}.d{j}*deg2m,PARAMS_bank{i}.z{j},'*b');
plot(PARAMS_banksym{i}.d{j}*deg2m,PARAMS_banksym{i}.z{j},'or');
%%
cs_length = 200; % meter
cs_delta = 1; % meter

i = 1; XYq = [640622.118316888,3495491.92378708];
% i = 3; XYq = [633852.118316888,3492831.92378708];

% [csd,csz,csx,csy,px,py] = GetCrossSections([X10(:),Y10(:),Z10(:)],SP{i},XYq,cs_length,cs_delta);
fDEM = griddedInterpolant(X',Y',Z');

i = 8;
PARAMS = [];

for j = 2 : length(Centerline{i})-1
    XYq = Centerline{i}(j,:);
[csd,csz,csx,csy,px,py] = GetCrossSectionsLinearInterp(fDEM,Centerline{i},XYq,cs_length,cs_delta);

% figure;
% plot(csd,csz,'k');

%%

Params = RetrieveCrossSectionParameters(csd,csz,[50 60]);

% Plot result
figure; hold on;
plot(csd,csz,'k');
tw1 = Params.TW/2;
x = [-tw1:tw1];
fz_p = @(TW,BW,H,x) (x < -BW/2 & x >= -TW/2).*(-2*H/(TW-BW)).*(x + BW/2)...
                + (x > BW/2 & x <= TW/2).*(2*H/(TW-BW)).*(x - BW/2);
plot(x+Params.CSC,fz_p(Params.TW,Params.BW,Params.H,x) + min(csz),'b','linewidth',1.5);


% store results
Params.x = XYq(1);
Params.y = XYq(2);
if isempty(PARAMS)
    PARAMS = Params;
else
    PARAMS(end+1) = Params;
end

end
PARAMS = struct2table(PARAMS);
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












