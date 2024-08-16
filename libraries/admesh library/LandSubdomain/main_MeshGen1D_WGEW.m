addpath(genpath('D:\academic\PR-PREEVENTS\codes_MATLAB\ADMESH\MeshGen1D'));
DEMFile = 'D:\academic\data\gis\walnut_gulch\dem\usgs_13arc_Merge_Clip.tif';
% DEMFile = 'D:\academic\data\gis\walnut_gulch\dem\usgs_1m_Merge_Clip_GCS_NAD83.tif';
BoundaryFile = 'D:\academic\data\gis\walnut_gulch\hydrography\WBD_AOI.shp';
NHDFile = 'D:\academic\data\gis\walnut_gulch\hydrography\NHDFlowline_AOI_2D.shp';

%%
Setting1D.h_min     = 30*km2deg(1e-3);
Setting1D.h_max     = 5e2*km2deg(1e-3);
Setting1D.h0        = 1*Setting1D.h_min;
Setting1D.K         = 20;
Setting1D.g         = 0.15;
Setting1D.SmoothingRMSE = 5*km2deg(1e-3); % set "0" to disable smoothing

%% ========================================================================
% Read data and smoothing
%==========================================================================
%--------------------------------------------------------------------------
% Extract Open-channels using TopoToolbox
%--------------------------------------------------------------------------
[FL0,FD,A] = ExtractOpenChannels(DEMFile);
k = 3e4;
FL1 = klargestconncomps(STREAMobj(FD,A > k));

%--------------------------------------------------------------------------
% Read watershed boundary
%--------------------------------------------------------------------------
BP = shaperead(BoundaryFile);
BP = [BP.X(:), BP.Y(:)];

%--------------------------------------------------------------------------
% Smoothing using TopoToolbox (CRS algorithm)
%--------------------------------------------------------------------------
k = 20;
SmoothFL = SmoothingCRS(FL1,BP,k);

[sx,sy] = STREAMobj2XY(FL1);
FL = NaNdlm2struct([sx,sy],'Boundary',BP);
FL2 = cellfun(@(x) [x; nan(1,2)],FL,'UniformOutput',0);
FL2 = vertcat(FL2{:});

figure; hold on; axis equal;
plot(BP(:,1),BP(:,2),'k');
plot(FL2(:,1),FL2(:,2));
% plot(sx,sy);

%
NHDFile = 'D:\academic\data\gis\walnut_gulch\hydrography\NHDFlowline_AOI_2D.shp';
SHP = shaperead(NHDFile);
SHP = struct2table(SHP);

SHP(SHP.Visibility < 500000,:) = [];

FL = [];
for i = 1 : size(SHP,1)
    x1 = SHP.X{i};
    y1 = SHP.Y{i};
    
    I = ~isnan(x1) & ~isnan(y1);
    I = I & inpolygon(x1,y1,BP(:,1),BP(:,2));

    x1 = x1(I);
    y1 = y1(I);
    
    FL{i} = [x1(:), y1(:)];
    
end

%% Compute projection functions of path
Points = FL;
% Points = MA.XY(52);

%----------------------------------------------------------------------
% Add fixed points for junctions
%----------------------------------------------------------------------
fixedPoints = [];
nfixedP = 0;
for j = 1 : length(Points)
    nfixedP = nfixedP + 1;
    fixedPoints(nfixedP,:) = [Points{j}(1,1),Points{j}(1,2)];
    nfixedP = nfixedP + 1;
    fixedPoints(nfixedP,:) = [Points{j}(end,1),Points{j}(end,2)];
end

clear PI;
k = 0;
for i = 1 : length(Points)
    x = Points{i}(:,1);
    y = Points{i}(:,2);

    PI1 = ComputePathCurvature([x,y],fixedPoints,Setting1D.SmoothingRMSE);
    if ~isempty(PI1)
        k = k + 1;
        PI(k) = PI1;
    end
end


figure; hold on; view(2);
for i = 1 : length(PI)
    plot(PI(i).x1(PI(i).p),PI(i).y1(PI(i).p),'k');
%     myPlot3(PI(i).sx(PI(i).p),PI(i).sy(PI(i).p),PI(i).k(PI(i).p),'linewidth',2);
%     myPlot3(PI(i).x1(PI(i).p),PI(i).y1(PI(i).p),deg2km(1./(20*PI(i).k(PI(i).p)))*1e3,'linewidth',2);
end
axis equal;

for i = 1 : length(PI)
    fprintf('%f %f\n',min(1./PI(i).k(PI(i).p))/Setting1D.h_min,max(1./PI(i).k(PI(i).p)));
end
%% Construct dummy constraints for ADMESH
ADMESH_folder = 'D:\academic\PR-PREEVENTS\codes_MATLAB\ADMESH\ADMESH_v3\';
Filename = [ADMESH_folder,'test_files\WGEW\GCS_NAD83\WGEW_GCS_NAD83_1m_linear'];
%----------------------------------------------------------------------
% Setup projection to use
%----------------------------------------------------------------------
for i = 1 : length(PI)
    PI(i).x = PI(i).x1;
    PI(i).y = PI(i).y1;
end

%----------------------------------------------------------------------
% Setup ADMESH input
%----------------------------------------------------------------------
PTS = [];

external_boundary = {BP};
for i = 1 : length(external_boundary)
    PTS.Poly(i).x = external_boundary{i}(:,1);
    PTS.Poly(i).y = external_boundary{i}(:,2);
end

for i = 1 : length(PI)
    p = PI(i).p;
    x1 = PI(i).x(p);
    y1 = PI(i).y(p);
    k1 = PI(i).k(p);
    
    PTS.Constraints(i).num = 18;
    PTS.Constraints(i).xy = [x1,y1];
    PTS.Constraints(i).type = 'line';
    PTS.Constraints(i).data = [];
    PTS.Constraints(i).Kappa = k1;
end

Settings = [];
Settings.Res = 'Low';
Settings.hmax = Setting1D.h_max;
Settings.hmin = Setting1D.h_min;
Settings.K.Status = 'on';
Settings.K.Value = Setting1D.K;

Settings.G.Status = 'on';
Settings.G.Value = Setting1D.g/sqrt(1);
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

%% Mesh generation 1D with updated mesh size
% Points = MA.XY(52);

%----------------------------------------------------------------------
% Add fixed points for junctions
%----------------------------------------------------------------------
fixedPoints = [];
nfixedP = 0;
for j = 1 : length(PI)
    nfixedP = nfixedP + 1;
    x = PI(j).x1(PI(j).p);
    y = PI(j).y1(PI(j).p);
    fixedPoints(nfixedP,:) = [x(1), y(1)];
    nfixedP = nfixedP + 1;
    fixedPoints(nfixedP,:) = [x(end), y(end)];
end
    
clear Mesh1D;
OverwriteFile = [Filename,'_dummyconstraint_updated_1d_hsize'];

for i = 1 : length(PI)
    x = PI(i).x1(PI(i).p);
    y = PI(i).y1(PI(i).p);

    id = ismember([x,y],fixedPoints,'rows');
    fixedPoints1D = PI(i).p(id);
    
    %----------------------------------------------------------------------
    % Generate 1D mesh
    %----------------------------------------------------------------------
    mesh1d = MeshGeneration1D(PI(i),fixedPoints1D,Setting1D,OverwriteFile);
    Mesh1D(i) = mesh1d;
end
% Mesh1D.Settings = Setting1D;

figure; hold on; axis equal;
% plot(x_b,y_b,'k');
for i = 1 : length(Mesh1D)
    x = PI(i).x1(PI(i).p);
    y = PI(i).y1(PI(i).p);
    plot(x,y,'k');
    plot(Mesh1D(i).X,Mesh1D(i).Y,'-*b');
end


for i = 1 : length(Mesh1D)
    fprintf('%f %f\n',min(diff(Mesh1D(i).p))/Setting1D.h_min,max(diff(Mesh1D(i).p))/Setting1D.h_max);
end
% for i = 1 : length(Points)
%     dx = diff(Mesh1D(i).X);
%     dy = diff(Mesh1D(i).Y);
%     d = sqrt(dx.^2 + dy.^2);
%     max(d)/Setting1D.h_max
% end

%% Write ADMESH input for single-line constraint (pass both 1D lines and 2D boundaries as constraints)
xl = min(BP(:,1)) - Setting1D.h_min;
xr = max(BP(:,1)) + Setting1D.h_min;
yl = min(BP(:,2)) - Setting1D.h_min;
yr = max(BP(:,2)) + Setting1D.h_min;

PTS = [];
external_boundary = {BP};

for i = 1 : length(external_boundary)
    PTS.Poly(i).x = external_boundary{i}(:,1);
    PTS.Poly(i).y = external_boundary{i}(:,2);
end

for i = 1 : length(PI)
    PTS.Constraints(i).num = 18;
    PTS.Constraints(i).xy = [Mesh1D(i).X,Mesh1D(i).Y];
    PTS.Constraints(i).type = 'line';
    PTS.Constraints(i).data = [];
    PTS.Constraints(i).Kappa = Mesh1D(i).K;
end


Settings = [];
Settings.Res = 'Low';
Settings.hmax = Setting1D.h_max;
Settings.hmin = Setting1D.h_min;
Settings.K.Status = 'on';
Settings.K.Value = Setting1D.K;

Settings.G.Status = 'on';
Settings.G.Value = Setting1D.g/sqrt(1);
Settings.View.Status = 'on';

Settings.R.Status = 'off';
Settings.B.Status = 'off';
Settings.T.Status = 'off';

Settings.DummyConstraint = 0;

% cd(ADMESH_folder);
ConstraintFilename = [Filename,'_single_constraint'];
save(ConstraintFilename,'PTS','Settings','Setting1D');




%% Retrieve bankline
% Read 10m-DEM file in matrix format
TIF = ReadTifFile(DEMFile);
X = TIF.X;
Y = TIF.Y;
Z = TIF.Z;

fDEM = griddedInterpolant(X',Y',Z');
fDEM = @(x,y) double(fDEM(x,y));

%%
clear BankLine CS;
for i = 1 : length(Mesh1D)
    cs_delta = 0.1*km2deg(1e-3);
    ChannelElem2DWidth = km2deg(100*1e-3);
    Centerline = [Mesh1D(i).X,Mesh1D(i).Y];
    [PARAMS_bank, CS_bankline] = RetrieveBankline(fDEM,Centerline,ChannelElem2DWidth,cs_delta,5,'degree');
    
    PARAMS_bank = table2struct(PARAMS_bank,'ToScalar',1);
    BankLine(i) = PARAMS_bank;
    
    CS_bankline = struct2table(CS_bankline);
    CS_bankline = table2struct(CS_bankline,'ToScalar',1);
    CS(i) = CS_bankline;
end

for i = 1 : length(Mesh1D)
    BankLine(i).X = [BankLine(i).Bxr; flip(BankLine(i).Bxl)];
    BankLine(i).Y = [BankLine(i).Byr; flip(BankLine(i).Byl)];
    BankLine(i).K = [Mesh1D(i).K; flip(Mesh1D(i).K)];
end


figure; hold on; axis equal;
% plot(x_b,y_b,'k');
for i = 1 : length(Mesh1D)
    x = PI(i).x1(PI(i).p);
    y = PI(i).y1(PI(i).p);
    plot(x,y,'k');
    plot(Mesh1D(i).X,Mesh1D(i).Y,'-*b');
    plot(BankLine(i).X,BankLine(i).Y,'k');
end


%% Debug
XYq = cursor_info.Position([1 2]);
Dist = inf;
for i = 1 : length(Mesh1D)
    Centerline = [Mesh1D(i).X,Mesh1D(i).Y];
    [j,dist] = knnsearch(Centerline,XYq);
    if dist < Dist
        Dist = dist;
        I = i;
        J = j;
    end
end

% j = 1794;
% j = 431;
figure; hold on;
% plot(csd,csz,'k');
plot(CS(I).d{J},CS(I).z{J},'k');
plot(BankLine(I).d{J},BankLine(I).z{J},'o');

%% Construct polygon from banklines
pgon = polyshape();
for i = 1 : length(Mesh1D)
    ipgon = polyshape(BankLine(i).X,BankLine(i).Y,'KeepCollinearPoints',true);
    pgon = union(pgon,ipgon,'KeepCollinearPoints',true);
end

% Remove holes: temporary solution
pgon = rmholes(pgon); 

% Remove duplicated vertices with tolerance 1 meter
V = pgon.Vertices;
while 1
    dV = diff(V);
    d = sqrt(sum(dV.^2,2));
    % I = find(d < km2deg(1e-3));
    I = find(d < Setting1D.h_min*0.5);
%     J = diff(I);
%     I = I(J > 1);
    if isempty(I)
        break;
    end
    V(I,:) = [];
end

% Vertices = uniquetol(pgon.Vertices,km2deg(0.1e-3),'byrows',true);

pgon = polyshape(V,'KeepCollinearPoints',true);


figure; hold on; axis equal;
% plot(x_b,y_b,'k');
for i = 1 : length(PI)
    x = PI(i).x1(PI(i).p);
    y = PI(i).y1(PI(i).p);
    plot(x,y,'k');
    plot(Mesh1D(i).X,Mesh1D(i).Y,'-*b');
end
plot(pgon);

%% Write ADMESH input (pass both 1D lines and 2D boundaries as constraints)
xl = min(BP(:,1)) - Setting1D.h_min;
xr = max(BP(:,1)) + Setting1D.h_min;
yl = min(BP(:,2)) - Setting1D.h_min;
yr = max(BP(:,2)) + Setting1D.h_min;

PTS = [];
external_boundary = {BP};

for i = 1 : length(external_boundary)
    PTS.Poly(i).x = external_boundary{i}(:,1);
    PTS.Poly(i).y = external_boundary{i}(:,2);
end

CenterXY = [vertcat(Mesh1D.X),vertcat(Mesh1D.Y)];
BankXY = [vertcat(BankLine.X),vertcat(BankLine.Y)];
BankK = vertcat(BankLine.K);
xy = pgon.Vertices;
xy(end+1,:) = xy(1,:);

id = knnsearch(BankXY,xy);
k = BankK(id);

PTS.Constraints.num = 18;
PTS.Constraints.type = 'line';
PTS.Constraints.data = [];
PTS.Constraints.xy = xy;
PTS.Constraints.Kappa = k;

Settings = [];
Settings.Res = 'Low';
Settings.hmax = Setting1D.h_max;
Settings.hmin = Setting1D.h_min;
Settings.K.Status = 'on';
Settings.K.Value = Setting1D.K;

Settings.G.Status = 'on';
Settings.G.Value = Setting1D.g/sqrt(1);
Settings.View.Status = 'on';

Settings.R.Status = 'off';
Settings.B.Status = 'off';
Settings.T.Status = 'off';

Settings.DummyConstraint = 0;

cd(ADMESH_folder);
ConstraintFilename = [Filename,'_bankline_constraint'];
save(ConstraintFilename,'PTS','Settings');







