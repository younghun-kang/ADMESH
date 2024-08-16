function PlotMesh(app,per)
% PlotMesh - Plots mesh
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% August 2013; Last revision: 12-August-2013

%--------------------------- BEGIN CODE ---------------------------------------

%------------------------------------------------------------------------------
% Get plot window handle
%------------------------------------------------------------------------------
pH = app.UIAxes; % Plot Window handle
disableDefaultInteractivity(pH); % appdesigner has an issue

%------------------------------------------------------------------------------
% If there is something currently in the plot window, delete it
%------------------------------------------------------------------------------
plotItems = get(pH, 'Children');

if ~isempty(plotItems); delete(plotItems); end

colorbar('delete') % Delete current colorbar

%------------------------------------------------------------------------------
% Tell ADMESH where to plot
%------------------------------------------------------------------------------
hold(pH,'on');

%------------------------------------------------------------------------------
% Set axes scaling
%------------------------------------------------------------------------------
nodes = app.MESH.Points(app.MESH.ConnectivityList,:);

xmin = min(nodes(:,1));
xmax = max(nodes(:,1));
ymin = min(nodes(:,2));
ymax = max(nodes(:,2));

% Compute image offset
offset = max((xmax-xmin),(ymax-ymin));

% Re-compute limits with offset
XMIN = xmin - offset*per;
XMAX = xmax + offset*per;
YMIN = ymin - offset*per;
YMAX = ymax + offset*per;

% Set plot ratios to auto mode
set(pH,'PlotBoxAspectRatioMode','auto','DataAspectRatioMode','auto')

% Set plot axis limits
axis(pH,[XMIN XMAX YMIN YMAX])

%------------------------------------------------------------------------------
% Plot mesh
%------------------------------------------------------------------------------
verts = app.MESH.Points;                        % vertices
faces = app.MESH.ConnectivityList;              % connectivity list

trisurf(faces,verts(:,1),verts(:,2),zeros(length(verts),1),...
    'FaceColor','none',...
    'EdgeColor',[0 .4 .8],...
    'Parent',pH,'Tag','Mesh');

% Set up right-click menu
% hcmenu = uicontextmenu;
% uimenu(hcmenu,'Label','Delete Element','Callback',@DeleteElement);
% set(pH,'tag','Mesh','uicontextmenu',hcmenu)

%------------------------------------------------------------------------------
% Plot boundaries
%------------------------------------------------------------------------------
nConstraints = numel(app.MESH.Constraints);

for k = 1:nConstraints

    if app.MESH.Constraints(k).num == -1 % Open
                
        Nodes = app.MESH.Constraints(k).nodeStr;

        h = line(verts(Nodes,1),verts(Nodes,2),'color','b','linewidth',3,'Parent',app.UIAxes);
        
        set(h,'tag','Mesh Constraint')
        
    elseif any(app.MESH.Constraints(k).num == [0 2 10 12 20 22 30]); % External Boundary
        
        Nodes = app.MESH.Constraints(k).nodeStr;

        h = line(verts(Nodes,1),verts(Nodes,2),'color','k','linewidth',3,'Parent',app.UIAxes);
        
        set(h,'tag','Mesh Constraint')
        
    elseif any(app.MESH.Constraints(k).num == [1 11 21]); % Internal Boundary
        
        Nodes = app.MESH.Constraints(k).nodeStr;

        h = line(verts(Nodes,1),verts(Nodes,2),'color',[0 .5 0],'linewidth',3,'Parent',app.UIAxes);
        
        set(h,'tag','Mesh Constraint')
        
    elseif any(app.MESH.Constraints(k).num == [3 13 23]); % External Constraints
                
        Nodes = app.MESH.Constraints(k).nodeStr;

        h = line(verts(Nodes,1),verts(Nodes,2),'color','r','linewidth',3,'Parent',app.UIAxes);
        
        set(h,'tag','Mesh Constraint')

    elseif any(app.MESH.Constraints(k).num == [18 17 19]); % Channel Constraints
                
        Nodes = app.MESH.Constraints(k).nodeStr;

        h = line(verts(Nodes,1),verts(Nodes,2),'color','r','linewidth',2,'Parent',app.UIAxes);
        
        set(h,'tag','Mesh Constraint')
        
    elseif any(app.MESH.Constraints(k).num == [4 5 24 25]); % Internal Constraints
                        
        % Plot first side
        Nodes = app.MESH.Constraints(k).nodeStr(:,1);

        h = line(verts(Nodes,1),verts(Nodes,2),'color','r','linewidth',3,'Parent',app.UIAxes);
        
        set(h,'tag','Mesh Constraint')
        
       % Plot second side
        Nodes = app.MESH.Constraints(k).nodeStr(:,2);
        
        h = line(verts(Nodes,1),verts(Nodes,2),'color','r','linewidth',3,'Parent',app.UIAxes);
        
        set(h,'tag','Mesh Constraint')

    end
    
end

% uistack(meshPatch,'top')

%------------------------------------------------------------------------------
% Show elevation if it is enabled
%------------------------------------------------------------------------------
%PlotElevation;

%------------------------------------------------------------------------------
% Set aspect ratio
%------------------------------------------------------------------------------
set(pH,'PlotBoxAspectRatio',[1 1 1],'DataAspectRatio',[1 1 1])
axis(pH,'tight')

%------------------------------------------------------------------------------
% Plot distance scale
%------------------------------------------------------------------------------
scalebar('hAxes',pH,'Location','southeast','Bold',1,'Unit','m')

%------------------------------------------------------------------------------
% Display mesh info
%------------------------------------------------------------------------------
DisplayMeshInfo(app);

%------------------------------------------------------------------------------
% Save axis limits
%------------------------------------------------------------------------------
app.xLimits = get(pH,'xlim');    % x-axis limits
app.yLimits = get(pH,'ylim');    % y-axis limits

end
