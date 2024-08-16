addpath(genpath('D:\academic\PR-PREEVENTS\codes_MATLAB\ADMESH\MeshGen1D'));
DEMFile = 'D:\academic\data\gis\middle_bosque_river\dem\usgs_3dep_13arc_Clip.tif';
BoundaryFile = 'D:\academic\data\gis\middle_bosque_river\hydrography\WBD_AOI.shp';
NHDFile = 'D:\academic\data\gis\middle_bosque_river\hydrography\NHDFlowline_AOI_2d.shp';
%%
Setting1D.h_min     = 200*km2deg(1e-3);
Setting1D.h_max     = 1000*km2deg(1e-3);
Setting1D.h0        = 1*Setting1D.h_min;
Setting1D.K         = 1;
Setting1D.g         = 0.15;
Setting1D.SmoothingRMSE = 10*km2deg(1e-3); % set "0" to disable smoothing

%% ========================================================================
% Read data and smoothing
%==========================================================================
%--------------------------------------------------------------------------
% Extract Open-channels using TopoToolbox
%--------------------------------------------------------------------------
[FL0,FD,A] = ExtractOpenChannels(DEMFile);
k = 1e5;
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


SHP = shaperead(NHDFile);
SHP = struct2table(SHP);
SHP(SHP.visibility < 500000,:) = [];
temp_x = cellfun(@(x) [x(:);nan],SHP.X,'UniformOutput',0);
temp_y = cellfun(@(x) [x(:);nan],SHP.Y,'UniformOutput',0);
temp_x = vertcat(temp_x{:});
temp_y = vertcat(temp_y{:});
plot(temp_x,temp_y);


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
for i = 1 : length(Points)
    x = Points{i}(:,1);
    y = Points{i}(:,2);

    PI(i) = ComputePathCurvature([x,y],fixedPoints,Setting1D.SmoothingRMSE);
end

figure; hold on;
for i = 1 : length(PI)
    plot(PI(i).x1(PI(i).p),PI(i).y1(PI(i).p),'k','linewidth',1.5);
    plot(PI(i).sx(PI(i).p),PI(i).sy(PI(i).p),'--b','linewidth',1.5);
%     myPlot3(PI(i).sx(PI(i).p),PI(i).sy(PI(i).p),PI(i).k(PI(i).p));
%     myPlot3(PI(i).sx(PI(i).p),PI(i).sy(PI(i).p),PI(i).k(PI(i).p),'LineWidth',2);
end
axis equal;
%% Construct pre-constraints for ADMESH
ADMESH_folder = 'D:\academic\PR-PREEVENTS\codes_MATLAB\ADMESH\ADMESH_v3\';
Filename = [ADMESH_folder,'test_files\MBRW\GCS_NAD83\MBRW_GCA_NAD83_hmin200hmax1000'];
%----------------------------------------------------------------------
% Setup projection to use
%----------------------------------------------------------------------
for i = 1 : length(Points)
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

for i = 1 : length(Points)
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
folder = fileparts(Filename);
cd(ADMESH_folder);
if ~exist(folder,'dir')
    mkdir(folder);
end
save(DummyConstFile,'PTS','Settings');


%% Mesh generation 1D with updated mesh size
Points = FL;
% Points = MA.XY(52);

load(OverwriteFile,'X','Y','h');
figure; hold on; daspect([1 1 1]);
temp = struct2table(PI);
temp = table2struct(temp,'ToScalar',1);
temp_x = cellfun(@(x,y) [x(y);nan],temp.x,temp.p,'UniformOutput',0);
temp_y = cellfun(@(x,y) [x(y);nan],temp.y,temp.p,'UniformOutput',0);
temp_x = vertcat(temp_x{:});
temp_y = vertcat(temp_y{:});
temp_h = interp2(X,Y,h,temp_x,temp_y);
myPlot3(temp_x,temp_y,temp_h,'linewidth',2);
caxis([Setting1D.h_min,Setting1D.h_max]);


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
    
clear Mesh1D;
OverwriteFile = [Filename,'_dummyconstraint_updated_1d_hsize'];
for i = 1 : length(Points)
    x = Points{i}(:,1);
    y = Points{i}(:,2);

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
for i = 1 : length(Points)
    plot(Points{i}(:,1),Points{i}(:,2),'k');
    plot(Mesh1D(i).X,Mesh1D(i).Y,'-*b');
end

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

id = 1 : length(Points);
for i = 1 : length(id)
    j = id(i);
    PTS.Constraints(j).num = 18;
    PTS.Constraints(j).xy = [Mesh1D(j).X,Mesh1D(j).Y];
    PTS.Constraints(j).type = 'line';
    PTS.Constraints(j).data = [];
    PTS.Constraints(j).Kappa = Mesh1D(j).K;
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
save(ConstraintFilename,'PTS','Settings');




%% Retrieve bankline
% Read 10m-DEM file in matrix format
TIF = ReadTifFile(DEMFile);
X = TIF.X;
Y = TIF.Y;
Z = TIF.Z;

fDEM = griddedInterpolant(X',Y',Z');


BankLine = [];
for i = 1 : length(Mesh1D.X)
    cs_delta = 0.1*km2deg(1e-3);
    ChannelElem2DWidth = 100*km2deg(1e-3);
    Centerline = [Mesh1D.X{i},Mesh1D.Y{i}];
    [PARAMS_bank, CS_bankline] = RetrieveBankline(fDEM,Centerline,ChannelElem2DWidth,cs_delta,2,'degree');
    
    BankLine.X{i} = [PARAMS_bank.Bxr; flip(PARAMS_bank.Bxl)];
    BankLine.Y{i} = [PARAMS_bank.Byr; flip(PARAMS_bank.Byl)];
    BankLine.K{i} = [Mesh1D.K{i}; flip(Mesh1D.K{i})];
end


%% Construct polygon from banklines
pgon = polyshape();
for i = 1 : length(Mesh1D.X)
    ipgon = polyshape(BankLine.X{i},BankLine.Y{i},'KeepCollinearPoints',true);
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
for i = 1 : length(Points)
    plot(Points{i}(:,1),Points{i}(:,2),'k');
    plot(Mesh1D.X{i},Mesh1D.Y{i},'-*b');
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

CenterXY = [vertcat(Mesh1D.X{:}),vertcat(Mesh1D.Y{:})];
BankXY = [vertcat(BankLine.X{:}),vertcat(BankLine.Y{:})];
BankK = vertcat(BankLine.K{:});
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






