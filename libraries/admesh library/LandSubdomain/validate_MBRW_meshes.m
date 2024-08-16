% addpath(genpath('D:\academic\PR-PREEVENTS\codes_MATLAB\MeshGeneration1D_v2'));
DEMFile = 'D:\academic\data\gis\middle_bosque_river\dem\3dep_13arc_Clip.tif';
BoundaryFile = 'D:\academic\data\gis\middle_bosque_river\hydrography\WBD_AOI.shp';

ADMESH_folder = 'D:\academic\PR-PREEVENTS\codes_MATLAB\ADMESH_yk\';
Filename = [ADMESH_folder,'test_files\MBRW\GCS_NAD83\MBRW_GCA_NAD83'];

SingleConstraintFilename = [Filename,'_single_constraint'];
BanklineConstraintFilename = [Filename,'_bankline_constraint'];
SmoothSingleConstraintFilename = [Filename,'_smooth_single_constraint'];
SmoothBanklineConstraintFilename = [Filename,'_smooth_bankline_constraint'];

MeshS = Read14File([SingleConstraintFilename,'.14']);
MeshB = Read14File([BanklineConstraintFilename,'.14']);
MeshSS = Read14File([SmoothSingleConstraintFilename,'.14']);
MeshSB = Read14File([SmoothBanklineConstraintFilename,'.14']);


%%
% load(SingleConstraintFilename);
load(SmoothSingleConstraintFilename);
mesh = MeshSS;

for i = 1 : length(PTS.Constraints)
Mesh1D.X{i} = PTS.Constraints(i).xy(:,1);
Mesh1D.Y{i} = PTS.Constraints(i).xy(:,2);
Mesh1D.K{i} = PTS.Constraints(i).Kappa;
end


%% Retrieve bankline
% Read 10m-DEM file in matrix format
TIF = ReadTifFile(DEMFile);
X = TIF.X;
Y = TIF.Y;
Z = TIF.Z;

fDEM = griddedInterpolant(X',Y',Z');


BankLine = [];
BanklineData = [];
CrossSectionData = [];
for i = 1 : length(Mesh1D.X)
    cs_delta = 0.1*km2deg(1e-3);
    ChannelElem2DWidth = 100*km2deg(1e-3);
    Centerline = [Mesh1D.X{i},Mesh1D.Y{i}];
%     Centerline = PTS.Constraints(i).xy;
    [BanklineData{i}, CrossSectionData{i}] = RetrieveBankline(fDEM,Centerline,ChannelElem2DWidth,cs_delta,2,'degree');
    BankLine.X{i} = [BanklineData{i}.Bxr; flip(BanklineData{i}.Bxl)];
    BankLine.Y{i} = [BanklineData{i}.Byr; flip(BanklineData{i}.Byl)];
%     BankLine.K{i} = [Mesh1D.K{i}; flip(Mesh1D.K{i})];
end

%% Construct polygon from banklines
pgon = polyshape();
for i = 1 : length(BankLine.X)
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
    I = find(d < Settings.hmin*0.5);
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
plot(pgon);
for i = 1 : length(Mesh1D.X)
    plot(Mesh1D.X{i},Mesh1D.Y{i},'-*b');
end


%%
i = 1;
PARAMS_bank = BanklineData{i};
CS = CrossSectionData{i};
CSParams{i} = RetrieveCrossSectionParameters(PARAMS_bank,CS);

%%
i = 1;
PARAMS_bank = BanklineData{i};
CS = CrossSectionData{i};
XYq = cursor_info.Position([1 2]);
j = knnsearch([PARAMS_bank.x,PARAMS_bank.y],XYq);
% j = 1794;
% j = 431;
figure; hold on;
% plot(csd,csz,'k');
plot(deg2km(CS(j).d)*1e3,CS(j).z,'k');
plot(deg2km(PARAMS_bank.d{j})*1e3,PARAMS_bank.z{j},'*b');




%%
mesh = MeshSB;

figure;
AX = axes;
hold on;
x = [CSParams{i}.Txl; flipud(CSParams{i}.Bxl)];
y = [CSParams{i}.Tyl; flipud(CSParams{i}.Byl)];
z = [CSParams{i}.BEL + CSParams{i}.H; flipud(CSParams{i}.BEL)];
patch(x,y,z,'edgecolor','interp');

x = [flipud(CSParams{i}.Bxl); CSParams{i}.Bxr];
y = [flipud(CSParams{i}.Byl); CSParams{i}.Byr];
z = [flipud(CSParams{i}.BEL); CSParams{i}.BEL];
patch(x,y,z,'edgecolor','interp');

x = [CSParams{i}.Bxr; flipud(CSParams{i}.Txr)];
y = [CSParams{i}.Byr; flipud(CSParams{i}.Tyr)];
z = [CSParams{i}.BEL; flipud(CSParams{i}.BEL+CSParams{i}.H)];
patch(x,y,z,'edgecolor','interp');
mycmap('dem',caxis); colorbar;
axis equal;

AX(2) = axes;
myPlot14File(MeshB,'field','mesh');

myLinkAxes(AX);

%%
figure;
AX_dem = axes;
AX_dem(2) = axes;

axes(AX_dem(2));
hold on; axis equal;
myPlot14File(mesh,'field','mesh');

% plot(PARAMS.Txr,PARAMS.Tyr,'k');
% plot(PARAMS.Txl,PARAMS.Tyl,'k');

% hold on; axis equal;
% for i = 1 : length(Mesh1D.X)
%     plot(Mesh1D.X{i},Mesh1D.Y{i},'b');
% end
% plot(pgon,'FaceColor','none','EdgeColor','k');


% x = [CSParams.Txr; flipud(CSParams.Txl)];
% y = [CSParams.Tyr; flipud(CSParams.Tyl)];
% plot(x,y,'k');

%%
% figure(1);
% AX_dem = AX_sym;
ax = axis;
cax = caxis(AX_dem(1));
[I,J] = find(X > ax(1) & X < ax(2) & Y > ax(3) & Y < ax(4));
I = unique(I);
J = unique(J);

surf(AX_dem(1),X(I,J),Y(I,J),Z(I,J),'EdgeColor','none','FaceColor','interp');
view(AX_dem(1),2);
axis(AX_dem(1),'equal');
colorbar(AX_dem(1));
if isequal(cax,[0 1])
    cax = zlim(AX_dem(1));
end
mycmap('dem',AX_dem(1),cax);

myLinkAxes(AX_dem);
% 









