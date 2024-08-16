function resetView(varargin)


%------------------------------------------------------------------------------
% Get GUI data
%------------------------------------------------------------------------------
app = varargin{1};

pH = app.UIAxes; % Plot Window handle

% Get axis settings
xLIM = app.xLimits; yLIM = app.yLimits;

set(pH,'xlim', xLIM)
set(pH,'ylim', yLIM)
drawnow

return

set(gui.Window,'renderer','opengl')

%------------------------------------------------------------------------------
% Check if there is anything currently in the plot window
%------------------------------------------------------------------------------
if isempty(get(pH, 'Children')) 
    return
end

set(gui.Window,'WindowButtonMotionFcn','')

% Get axis settings
xLIM = gui.xLimits; yLIM = gui.yLimits;

% Get current axis limits
pHxLIM = get(pH,'xlim'); pHyLIM = get(pH,'ylim');

% Bounding box
bbox = [pHxLIM pHyLIM];

% Determine if step is needed
dt = 20;

dxmin = linspace(bbox(1),xLIM(1),dt);
dxmax = linspace(bbox(2),xLIM(2),dt);
dymin = linspace(bbox(3),yLIM(1),dt);
dymax = linspace(bbox(4),yLIM(2),dt);

for i = 1:dt
    
    set(pH,'xlim', [dxmin(i) dxmax(i)])
    set(pH,'ylim', [dymin(i) dymax(i)])
    drawnow
    
end

% Check for google earth image
if ~isempty(findobj(get(pH,'children'),'tag','googlemap'))
    GetGoogleEarthImage(guiFig)
end

set(gui.Window,'renderer','zbuffer')

set(gui.Window,'WindowButtonMotionFcn',@CoordDisplay)

end