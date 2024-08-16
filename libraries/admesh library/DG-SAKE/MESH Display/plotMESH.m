%==========================================================================
% Plot mesh
%
% Author:
% Dustin West
% Assistant Engineer | Hazen and Sawyer               
% dwest@hazenandsawyer.com | hazenandsawyer.com
%==========================================================================

%--------------------------------------------------------------------------
% Open Figure
%--------------------------------------------------------------------------
fig = figure('color','w'); hold on; 

cax = get(gca,'colororder');

MESH2plot = MESH;
MESH2plot.Points(:,3) = -MESH2plot.Points(:,3);
% Get coordinates
x = MESH2plot.Points(:,1);
y = MESH2plot.Points(:,2);
z = MESH2plot.Points(:,3);

%--------------------------------------------------------------------------
% Compute colormap index for mesh
%--------------------------------------------------------------------------
zmin = min(z); % Compute min elevation
zmax = max(z); % Compute max elevation

% Normalize and multiply by 256 (number of rows in colormap)
C = round(256*(z - zmin) / ( zmax - zmin ));
C(C>256) = 256;  C(C <= 0) = 1; % Cap

% Get landcover colormap
cmap = landcolor(256);

%--------------------------------------------------------------------------
% Get MESH edges
%--------------------------------------------------------------------------
E = unique(...
    sort(...
    [MESH2plot.ConnectivityList(:,[1 2]); ...
    MESH2plot.ConnectivityList(:,[2 3]); ...
    MESH2plot.ConnectivityList(:,[3 1])],2),'rows');

%--------------------------------------------------------------------------
% Get Channel Routing Edges
%--------------------------------------------------------------------------
clear RE;
RE = cell(0);
try
catch
    RE = cell(0);
end
%--------------------------------------------------------------------------
% Get Channel Edges
%--------------------------------------------------------------------------
CE = cell(length(MESH2plot.Constraints),1);
for k = 1:length(MESH2plot.Constraints)
    if MESH2plot.Constraints(k).num == 18 || MESH2plot.Constraints(k).num == 6
        % Get coordinates
        if isfield(MESH2plot.Constraints,'nodes')
            if any(MESH2plot.Constraints(k).data(:,3) == 0)
                RE{end+1} = [MESH2plot.Constraints(k).nodes(1:end-1) MESH2plot.Constraints(k).nodes(2:end)];
            else
        CE{k} = [MESH2plot.Constraints(k).nodes(1:end-1) MESH2plot.Constraints(k).nodes(2:end)];
            end
        
        elseif isfield(MESH2plot.Constraints,'nodeStr')
            if any(MESH2plot.Constraints(k).data(:,3) == 0)
                RE{end+1} = [MESH2plot.Constraints(k).nodeStr(1:end-1) MESH2plot.Constraints(k).nodeStr(2:end)];
            else
                CE{k} = [MESH2plot.Constraints(k).nodeStr(1:end-1) MESH2plot.Constraints(k).nodeStr(2:end)];
            end
        end
    end
end
try
    RE = sort(vertcat(RE{:}),2);
catch
    RE = [];
end

CE = sort(vertcat(CE{:}),2);

%--------------------------------------------------------------------------
% Remove identified channel edges from routing edges
%--------------------------------------------------------------------------
if ~isempty(RE)
end
%--------------------------------------------------------------------------
% Remove all channel edges from E
%--------------------------------------------------------------------------
if ~isempty(RE); E = setdiff(E,RE,'rows'); end
if ~isempty(CE); E = setdiff(E,CE,'rows'); end

%--------------------------------------------------------------------------
% Plot MESH surface
%--------------------------------------------------------------------------
% patch(...
%     'faces',MESH.ConnectivityList,...
%     'vertices',MESH.Points,...
%     'facevertexcdata',cmap(C,:),...
%     'facecolor','interp',...
%     'edgecolor','none');

trisurf(MESH2plot.ConnectivityList,x,y,z,'facecolor','interp','edgecolor','none')

colors = get(gca,'colororder');

%--------------------------------------------------------------------------
% Plot element edges
%--------------------------------------------------------------------------
c = 60;
% line(x(E'),y(E'),z(E'),'color',[c c c]./255)

x1 = x(E'); x1(3,:) = nan; x1 = x1(:);
y1 = y(E'); y1(3,:) = nan; y1 = y1(:);
z1 = z(E'); z1(3,:) = nan; z1 = z1(:);
plot3(x1,y1,z1,'color',[c c c]./255);

%--------------------------------------------------------------------------
% Plot  channel edges
%--------------------------------------------------------------------------
% h1 = line(x(CE'),y(CE'),z(CE'),'color','b','linewidth',2);

x1 = x(CE'); x1(3,:) = nan; x1 = x1(:);
y1 = y(CE'); y1(3,:) = nan; y1 = y1(:);
z1 = z(CE'); z1(3,:) = nan; z1 = z1(:);
h1 = plot3(x1,y1,z1,'color','b','linewidth',2);

%--------------------------------------------------------------------------
% Plot  routing edges
%--------------------------------------------------------------------------
% h2 = line(x(RE'),y(RE'),z(RE'),'color','r','linewidth',2);

x1 = x(RE'); x1(3,:) = nan; x1 = x1(:);
y1 = y(RE'); y1(3,:) = nan; y1 = y1(:);
z1 = z(RE'); z1(3,:) = nan; z1 = z1(:);
h2 = plot3(x1,y1,z1,'color','r','linewidth',2);

%--------------------------------------------------------------------------
% Set axis properties
%--------------------------------------------------------------------------
colormap(cmap)
daspect([1 1 1/13])
minxyz = min(MESH2plot.Points,[],1);
maxxyz = max(MESH2plot.Points,[],1);
axis([minxyz(1) maxxyz(1) minxyz(2) maxxyz(2) minxyz(3) maxxyz(3)+50])
 
cbh         = colorbar('location','eastoutside');
cbh.YTick = linspace(0,.95,10);
set(cbh,'yticklabel',arrayfun(@(x) num2str(x,'%.0f'),round(linspace(minxyz(3), maxxyz(3),10)),'uni',false));
set(get(cbh,'Title'),'String','m')
xlabel('X (m)'), ylabel('Y (m)'), zlabel('Z (m)')
set(gca,'Xtick',linspace(minxyz(1),maxxyz(1),9))
set(gca,'Ytick',linspace(minxyz(2),maxxyz(2),5))
set(gca,'fontsize',18,'fontname','Times')
box on
grid on
view(3)
ax          = gca;                      % Axes handle
axpos       = ax.Position;              % Axes position
cpos        = cbh.Position;             % cbar position
cpos(3)     = 0.5*cpos(1,3);              % Thin bar
%cpos(4)     = 0.6*cpos(1,4);              % Thin bar
cbh.Position= cpos;                     % Set
ax.Position = axpos;                    % Re-set position
legend([h1(1) h2(1)],{'1D Channel Network','1D Overland Flow Routing'},...
    'location','southeast')

clear MESH2plot;
