function PlotEdgeStructure(varargin)
% PlotEdgeStructure - Plots the domain defined in PTS
%
% Syntax:  PlotEdgeStructure(PTS, IBtype, per, guiFig)
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
% Get GUI data
%------------------------------------------------------------------------------
if nargin == 2
    [app,per]   = deal(varargin{:});
%     [PTS,sb,pH] = deal(gui.PTS,[],[]);
    PTS = app.PTS;
    pH = app.UIAxes;
else
    [PTS,app,sb,pH,per] = deal(varargin{:}); % NOAA Coastline Input
end

%------------------------------------------------------------------------------
% If there is something currently in the plot window, delete it
%------------------------------------------------------------------------------

% Get all graphics object handles
plotItems = get(pH, 'Children'); 

% Delete all current graphics object handles
if ~isempty(plotItems); delete(plotItems); colorbar('delete'); end

% Remove mesh legend if it exist
% delete(findall(gui.Window,'Tag','Mesh Info')) 

%------------------------------------------------------------------------------
% Compute bounding box & set axis limits
%------------------------------------------------------------------------------

% Compute limits
xmin = min(vertcat(PTS.Poly(:).x));
xmax = max(vertcat(PTS.Poly(:).x));
ymin = min(vertcat(PTS.Poly(:).y));
ymax = max(vertcat(PTS.Poly(:).y));

% Compute image offset
offset = max((xmax-xmin),(ymax-ymin));

% Re-compute limits with offset
XMIN = xmin - offset*per;
XMAX = xmax + offset*per;
YMIN = ymin - offset*per;
YMAX = ymax + offset*per;

% Set plot ratios to auto mode
% set(pH,'PlotBoxAspectRatioMode','auto','DataAspectRatioMode','auto')

% Set plot axis limits
axis(pH,[XMIN XMAX YMIN YMAX])

msg = 'Plotting coastline...';
msg = 'Displaying final mesh...';
uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');
%------------------------------------------------------------------------------
% Initialize linewidth constant & graphics handle for user data
%------------------------------------------------------------------------------
LineWidth   = 1;
gH          = zeros(numel(PTS.Poly),1);

%------------------------------------------------------------------------------
% Tell ADMESH where to plot
%------------------------------------------------------------------------------
% axes(pH); 
hold(pH,'on');

%------------------------------------------------------------------------------
% Plot external boundary
%------------------------------------------------------------------------------
gH(1) = line(PTS.Poly(1).x , PTS.Poly(1).y,'Parent',pH);

set(gH(1),...
    'color','k',...
    'LineWidth',LineWidth,...
    'tag','Edge Structure',...
    'UserData',{'External Boundary',1})

%------------------------------------------------------------------------------
% Plot internal boundaries
%------------------------------------------------------------------------------
for k = 2:length(PTS.Poly)
    
    gH(k) = line(PTS.Poly(k).x , PTS.Poly(k).y,'Parent',pH);
    
    set(gH(k),...
        'color',[0 .5 0],...
        'LineWidth',LineWidth,...
        'tag','Edge Structure',...
        'UserData',{'Internal Boundary',k})
    
end

%------------------------------------------------------------------------------
% Define a context menu
%------------------------------------------------------------------------------

% Define the context menu items and install their callbacks for external boundary
hcmenu = uicontextmenu(app.UIFigure); % Create user interface context menu
uimenu(hcmenu,'Label','Assign Open Ocean Boundary Condition','MenuSelectedFcn',@(menu,action)AssignBC(app,menu,action));
% m1 = uimenu(hcmenu,'Text','Assign Open Ocean Boundary Condition','MenuSelectedFcn',@(mn) AssignBC(app,mn));
set(gH(1),'ContextMenu',hcmenu)

% Define the context menu items and install their callbacks for internal boundaries
hcmenu = uicontextmenu(app.UIFigure); % Create user interface context menu
%uimenu(hcmenu,'Label','Assign BC'       ,'Callback',@AssignBC);
uimenu(hcmenu,'Label','Delete Polygon'  ,'Callback',@DeletePolygon);

for k = 2:length(gH)
    set(gH(k),'uicontextmenu',hcmenu)
end

%------------------------------------------------------------------------------
% Plot constraints
%------------------------------------------------------------------------------
PlotConstraints(app,pH);

%------------------------------------------------------------------------------
% Set aspect ratio
%------------------------------------------------------------------------------
set(pH,'PlotBoxAspectRatio',[1 1 1],'DataAspectRatio',[1 1 1])
axis(pH,'tight');

%------------------------------------------------------------------------------
% Plot distance scale
%------------------------------------------------------------------------------
h = uiaxes(app.UIFigure);
scalebar('hAxes',app.UIAxes,'Location','southeast','Bold',1,'Unit','m');

drawnow; pause(.005)

%------------------------------------------------------------------------------
% Save limits
%------------------------------------------------------------------------------
app.xLimits = get(pH,'xlim');    % x-axis limits
app.yLimits = get(pH,'ylim');    % y-axis limits

end