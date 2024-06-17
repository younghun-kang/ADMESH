function status = ClearWindow(app)
% ClearAllCallback - GUI Callback that clears current mesh
% file.
%
% Syntax:  ClearMeshCallback(varargin)
%
% Inputs:
%    guiFig - handle that identifies the figure
%
% Outputs:
%    non
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% August 2013; Last revision: 08-August-2013

%------------------------- BEGIN CODE ------------------------------------
status = 1;

%--------------------------------------------------------------------------
% Get GUI data
%--------------------------------------------------------------------------
pH  = app.UIAxes; % Plot Window handle

obj1 = findobj(get(pH,'children'),'tag','Edge Structure');
obj2 = findobj(get(pH,'children'),'tag','Mesh');

if isempty(obj1) && isempty(obj2)
   return 
end

%--------------------------------------------------------------------------
% Ask user if they are sure they want to delete the mesh
%--------------------------------------------------------------------------
msg = 'Are you sure you want to clear?';
choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
    'Options',{'Clear','Cancel'},'DefaultOption',2,'Icon','Warning');

drawnow; pause(0.05);  % this innocent line prevents the Matlab hang

if strcmp(choice,'Cancel')
    status = 0;
    return
end

% Turn off colormap
% SetContourStatus(app.Window,nan,'off'); drawnow

%--------------------------------------------------------------------------
% Clear mesh
%--------------------------------------------------------------------------
if ~isempty(app.PTS) && ~isempty(obj2)
    
    plotItems = get(pH, 'Children');
    
    if ~isempty(plotItems); delete(plotItems); end
    
    colorbar('delete') % Delete current colorbar
    
    % Remove current legend if it exist
    app.ResultsBox.Value = '';
    
    app.MESH               = []; % Mesh data
    
    % guidata(fig,gui)
    
    PlotEdgeStructure(app,.1);
    
    return
end

%--------------------------------------------------------------------------
% Clear everything
%--------------------------------------------------------------------------
plotItems = get(pH, 'Children');

if ~isempty(plotItems); delete(plotItems); end

colorbar('delete') % Delete current colorbar

% Remove current legend if it exist
app.ResultsBox.Value = '';

SetContourStatus(app,'Off');

%--------------------------------------------------------------------------
% Clear app data
%--------------------------------------------------------------------------
app.PTS                = [];   % Edge Structure
app.xyzFun             = [];   % Elevation interpolant
app.MESH               = []; % Mesh data
app.MinEQ              = 0.3;  % Minimum element quality
app.xLimits            = [0 1];
app.yLimits            = [0 1];
app.per                = 0;    % Percent offset
app.ElevationDataFilename = [];

drawnow;

end