DEMFile = 'D:\academic\data\gis\walnut_gulch\dem\usgs_1m_Merge_Clip_GCS_NAD83.tif';
DEMFile = 'D:\academic\data\gis\walnut_gulch\dem\usgs_13arc_Merge_Clip.tif';

% DEMFile = 'D:\academic\data\gis\middle_bosque_river\dem\3dep_13arc_Clip_ProjectRast.tif';

%% Read 10m-DEM file in matrix format
TIF = ReadTifFile(DEMFile);
X = TIF.X;
Y = TIF.Y;
Z = TIF.Z;

%%
load('D:\academic\PR-PREEVENTS\codes_MATLAB\ADMESH\ADMESH_v3\test_files\WGEW\GCS_NAD83\WGEW_GCS_NAD83_1m_linear_single_constraint');

Mesh1D = [];
for i = 1 : length(PTS.Constraints)
    Mesh1D.X{i} = PTS.Constraints(i).xy(:,1);
    Mesh1D.Y{i} = PTS.Constraints(i).xy(:,2);
end
%%
fDEM = griddedInterpolant(X',Y',Z');
m2deg = km2deg(1e-3);

i = 1;
cs_delta = 0.1*km2deg(1e-3);
ChannelElem2DWidth = 100*m2deg;
Centerline = [Mesh1D.X{i},Mesh1D.Y{i}];
[PARAMS_bank, CS_bankline] = RetrieveBankline(fDEM,Centerline,ChannelElem2DWidth,cs_delta,2,'degree');
[PARAMS_banksym, CS_bankline] = RetrieveBankline(fDEM,Centerline,ChannelElem2DWidth,cs_delta,5,'degree');
[PARAMS_bank6, CS_bankline] = RetrieveBankline(fDEM,Centerline,ChannelElem2DWidth,cs_delta,6,'degree');

figure;
hold on;
plot(Centerline(:,1),Centerline(:,2),'k');
plot(PARAMS_bank.Bxr,PARAMS_bank.Byr,'b');
plot(PARAMS_bank.Bxl,PARAMS_bank.Byl,'b');
plot(PARAMS_banksym.Bxr,PARAMS_banksym.Byr,'r');
plot(PARAMS_banksym.Bxl,PARAMS_banksym.Byl,'r');
plot(PARAMS_bank6.Bxr,PARAMS_bank6.Byr,'color',mycolors('dgn'));
plot(PARAMS_bank6.Bxl,PARAMS_bank6.Byl,'color',mycolors('dgn'));
axis equal;

%% Debug
XYq = cursor_info.Position([1 2]);
j = knnsearch([PARAMS_bank.x,PARAMS_bank.y],XYq);
% j = 1794;
% j = 431;
figure; hold on;
% plot(csd,csz,'k');
plot(CS_bankline(j).d,CS_bankline(j).z,'k');
plot(PARAMS_bank.d{j},PARAMS_bank.z{j},'*b');
plot(PARAMS_banksym.d{j},PARAMS_banksym.z{j},'or');
plot(PARAMS_bank6.d{j},PARAMS_bank6.z{j},'o','color',mycolors('dgn'));

% figure;
% plot(CS_bankline(j).x,CS_bankline(j).y,'k');

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












