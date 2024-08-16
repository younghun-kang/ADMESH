%% Read Mesh
MBRW = Read14File('D:\academic\PR-PREEVENTS\codes_MATLAB\ADMESH_yk\ADMESH version 2.0_channel\Middle_Bosque_river_watershed_try3.14');


%% Get open channel boundaries

k = 0;
x = MBRW.TR.Points(:,1);
y = MBRW.TR.Points(:,2);
TRIx = x(MBRW.TR.ConnectivityList);
TRIy = y(MBRW.TR.ConnectivityList);

AvgElemSize = [];
for i = 1 : length(MBRW.BoundarySegments)
    if strcmpi(MBRW.BoundarySegments(i).Type,'NormalFluxType18')
        k = k + 1;
        SP{k} = [x(MBRW.BoundarySegments(i).Nodes), y(MBRW.BoundarySegments(i).Nodes)];
        
        ElemSize = zeros(length(MBRW.BoundarySegments(i).Nodes)-1,1);
        for j = 1 : length(MBRW.BoundarySegments(i).Nodes)-1
            id = ismember(TRIx,SP{k}([j j+1],1)) & ismember(TRIy,SP{k}([j j+1],2));
            id = find(sum(id,2) == 2);
            if length(id) ~= 2
                continue;
            end
            
            d = 0;
            for jj = 1 : length(id)
            iTRIx = TRIx(id(jj),:);
            iTRIy = TRIy(id(jj),:);
            
            jd = find(~ismember(TRIx(id(jj),:),SP{k}([j j+1],1)) | ~ismember(TRIy(id(jj),:),SP{k}([j j+1],2)));
            
            if jd == 1
                i1 = 1; i2 = 2; i3 = 3;
            elseif jd == 2
                i1 = 2; i2 = 1; i3 = 3;
            elseif jd == 3
                i1 = 3; i2 = 1; i3 = 2;
            else
                
            end
            
            a = [iTRIx(i2) - iTRIx(i3), iTRIy(i2) - iTRIy(i3)];
            b = [iTRIx(i1) - iTRIx(i3), iTRIy(i1) - iTRIy(i3)];
            d1 = dot(a,b)/norm(a);
            d1 = sqrt(norm(b).^2 - d1.^2);
            
            d = d + d1;
                
                
            end
            ElemSize(j) = d /2;
        end
        AvgElemSize{k} = [ElemSize(1); (ElemSize(1:end-1)+ElemSize(2:end))/2; ElemSize(end)];
    end
end

% figure; hold on;
% for i = 1 : length(SP)
%     myPlot3(SP{i}(:,1),SP{i}(:,2),'b','linewidth',1);
% end
% axis equal;

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
path = 'D:\academic\data\gis\archive\Bosque River Watershed';
file = [path,'\ned10m_MBRW.tif'];
TIF = ReadTifFile(file);
X10 = TIF.X;
Y10 = TIF.Y;
Z10 = TIF.Z;
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
% figure;
% surf(X10,Y10,Z10,'EdgeColor','none','FaceColor','interp');
% view(2); axis equal tight; mycmap('dem'); colorbar;

%%
fDEM = griddedInterpolant(X10',Y10',Z10');

i = 1;
cs_delta = 1;
ChannelElem2DWidth = 2*AvgElemSize{i};
[PARAMS, CS] = RetrieveCrossSectionParameters_v2(fDEM,SP{i},ChannelElem2DWidth,cs_delta);

% plot 2d
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

AXmesh(2) = axes;
myPlot14File(MBRW,'field','mesh');

myLinkAxes(AXmesh);

%% Debug
XYq = cursor_info.Position([1 2]);
j = knnsearch([PARAMS.x,PARAMS.y],XYq);
% j = 1794;
% j = 431;
figure; hold on;
% plot(csd,csz,'k');
tw1 = PARAMS.TW(j)/2;
bw1 = PARAMS.BW(j)/2;
x = [linspace(-tw1,-bw1),linspace(bw1,tw1)];
fz_p = @(TW,BW,H,x) (x < -BW/2 & x >= -TW/2).*(-2*H/(TW-BW)).*(x + BW/2)...
                + (x > BW/2 & x <= TW/2).*(2*H/(TW-BW)).*(x - BW/2);
plot(x,fz_p(PARAMS.TW(j),PARAMS.BW(j),PARAMS.H(j),x) + PARAMS.BEL(j),'b','linewidth',1.5);
plot(CS(j).d,CS(j).z,'k');

%%
fDEM = griddedInterpolant(X10',Y10',Z10');

i = 1;
cs_delta = 0.1;
ChannelElem2DWidth = 2*AvgElemSize{i};
ChannelElem2DWidth = 200;
[PARAMS_bank, CS_bankline] = RetrieveBankline(fDEM,SP{i},ChannelElem2DWidth,cs_delta);
[PARAMS_banksym, CS_bankline] = RetrieveBanklineSym(fDEM,SP{i},ChannelElem2DWidth,cs_delta);

figure;
AX_bankline = axes;
myPlot14File(MBRW,'field','mesh');

AX_bankline(2) = axes;
hold on;
plot(PARAMS_bank.Bxr,PARAMS_bank.Byr,'r');
plot(PARAMS_bank.Bxl,PARAMS_bank.Byl,'r');

myLinkAxes(AX_bankline);

%% Debug
XYq = cursor_info.Position([1 2]);
j = knnsearch([PARAMS_bank.x,PARAMS_bank.y],XYq);
% j = 1794;
% j = 431;
figure; hold on;
% plot(csd,csz,'k');
plot(CS_bankline(j).d,CS_bankline(j).z,'k');
plot(PARAMS_bank.d{j},PARAMS_bank.z{j},'o');
%%
cs_length = 200; % meter
cs_delta = 1; % meter

i = 1; XYq = [640622.118316888,3495491.92378708];
% i = 3; XYq = [633852.118316888,3492831.92378708];

% [csd,csz,csx,csy,px,py] = GetCrossSections([X10(:),Y10(:),Z10(:)],SP{i},XYq,cs_length,cs_delta);
fDEM = griddedInterpolant(X10',Y10',Z10');

i = 8;
PARAMS = [];

for j = 2 : length(SP{i})-1
    XYq = SP{i}(j,:);
[csd,csz,csx,csy,px,py] = GetCrossSectionsLinearInterp(fDEM,SP{i},XYq,cs_length,cs_delta);

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












