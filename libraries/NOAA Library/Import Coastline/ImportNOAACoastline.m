function ImportNOAACoastline(varargin)
% ImportCoastline -Uses the NOAA online coastline extractor to create a domain 
%
% Syntax:  ImportCoastline(guiFig,~,~)
%
% Inputs:  
%    guiFig - handle that identifies the figure
%
% Outputs:
%    Function produces a .mat file for ADmesh
%
% Other m-files required: none
% Subfunctions: 
% MAT-files required: none

% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% August 2013; Last revision: 08-August-2013

%------------- BEGIN CODE --------------

%--------------------------------------------------------------------------
% Get GUI data
%--------------------------------------------------------------------------
fig = varargin{1}; gui = % guidata(fig);

% Turn off colormap
SetContourStatus(gui.Window,nan,'off');

%--------------------------------------------------------------------------
% Ask user for filename
%--------------------------------------------------------------------------
[filename, pathname] = uigetfile( '*.txt*', 'Select a file');

% If user cancels
if ~filename; return; end

%--------------------------------------------------------------------------
% Open File ; Read in ; Check for consistency
%--------------------------------------------------------------------------

% Open file
fid = fopen([pathname filename],'r');

% Read in (lon,lat) data
try
    Coastline_Data = cell2mat(textscan(fid,'%f %f %*[^\n]'));
catch
    errordlg('There was an error readin in this file.','Error');
    return
end

% Close file
fclose(fid);

% Check data
if size(Coastline_Data,1) < 3
    errordlg('Not enough data points to create a domain.','Error');
    return
end

%--------------------------------------------------------------------------
% Reset variables
%--------------------------------------------------------------------------
gui.FilePath    = [];
gui.PTS         = [];
gui.xyzFun      = [];
gui.MESH        = [];

% update GUI workspce
% guidata(fig,gui);

%--------------------------------------------------------------------------
% If there is something currently in the plot window, delete it
%--------------------------------------------------------------------------
pH = gui.ViewAxes;

% Get all graphics object handles
plotItems = get(pH, 'Children'); 

% Delete all current graphics object handles
if ~isempty(plotItems); delete(plotItems); colorbar('delete'); end

% Remove mesh legend if it exist
delete(findall(gui.Window,'Tag','Mesh Info')) 

% Open NOAA coastline extractor website
%web http://www.ngdc.noaa.gov/mgg_coastline/mapit.jsp -browser

%--------------------------------------------------------------------------
% Build PTS Data Structure
%--------------------------------------------------------------------------
[PTS] = BuildEdgeStructure(Coastline_Data,pH,gui);

if isempty(PTS);
    
    % gui.sb.setText('Ready')

    set(pH ,'XLim',[0 1],'YLim',[0 1])
    drawnow
    
    return
end

% Convert coordinates to cartesian
% gui.sb.setText('Converting to cartesian coordinates...')

%--------------------------------------------------------------------------
% Add Edge Structure Attributes
%--------------------------------------------------------------------------
while 1
    
    % gui.sb.setText('Draw open ocean boundary...')
    
    %----------------------------------------------------------------------
    % Set axis prpoerties
    %----------------------------------------------------------------------
    PlotEdgeStructure(PTS,gui,gui.sb,pH,.6)

    zoomH = findobj('Tag','Zoom Tag');
    panH  = findobj('Tag','Pan Tag');
    
    % Set all other toolbar tools off
    set([zoomH panH],'State','off')
    
    % Have user select open ocean boundary
    tPTS = DrawOpenOceanBoundary(PTS,gui,pH,gui.sb);
    
    % Plot new PTS
    PlotEdgeStructure(tPTS,gui,gui.sb,pH,.1)

    % gui.sb.setText('Continue?')
    
    % Ask the user if they like their selection
    msg = [...
        'Select Continue to keep the coordinates of the open ocean '...
        'boundary you entered. Select Re-Try to enter different '...
        'coordinates.'];
    choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'Continue','Re-Try','Quit'},'DefaultOption',1,'Icon','Warning');
    drawnow; pause(0.05);  % this innocent line prevents the Matlab hang
    
    switch choice
        
        case 'Continue'
            
            PTS = tPTS; clear tPTS
            break
            
        case 'Re-Try'
            
            clear tPTS 
            
            % Clear window
            plotItems = get(gui.PlotWindow, 'Children');
            
            if ~isempty(plotItems); delete(plotItems); end
            
            PlotEdgeStructure(PTS,gui,gui.sb,pH,.6)
 
        case 'Quit'
            
            clear tPTS
            
            % Clear window
            plotItems = get(gui.PlotWindow, 'Children');
            
            if ~isempty(plotItems); delete(plotItems); end
            
            % gui.sb.setText('Ready')
            
            set(pH ,...
                'XLim',[0 1],...
                'YLim',[0 1])
            drawnow
            
            % gui.sb.setText('Ready')

            return
            
    end
    
    
end

%--------------------------------------------------------------------------
% Check for islands outside main domain
%--------------------------------------------------------------------------
% gui.sb.setText('Checking domain for exterior islands...')
gui.sb.ProgressBar.setVisible(true)
nSegments = numel(PTS.Poly);
set(gui.sb.ProgressBar, 'Minimum',1, 'Maximum',nSegments, 'Value',1)

ir = false(1,nSegments);

for i = 2:nSegments
    
    set(gui.sb.ProgressBar,'Value',i)
        
    IN = PointInPolygon(PTS.Poly(i).x, PTS.Poly(i).y,PTS.Poly(1).x, PTS.Poly(1).y);
    
    if any(~IN); ir(i) = true; end 
    
end

gui.sb.ProgressBar.setVisible(false) 

PTS.Poly(ir) = []; % Remove exterior islands

% gui.sb.setText('Save File...')

% Save as a .mat file
[filename,pathname] = uiputfile('*.mat','Save File As');

if filename ~=0
    
    xyzFun = []; 
        
    gui.FilePath = [pathname filename];% Save filename/pathname
    
    % gui.sb.setText('Saving File...')
    
    save([pathname, filename], 'PTS', 'xyzFun','-v7.3');
    
    gui.FilePath = [pathname filename];
    
    % Update GUI workspace
    gui.PTS     = PTS;
    gui.xyzFun  = xyzFun;
    
    % gui.sb.setText('Plotting in cartesian coordinates')

    % Plot Boundary
    PlotEdgeStructure(gui,.1);
    
    % gui.sb.setText( 'Ready')
    
else
    
    % Clear window
    plotItems = get(gui.PlotWindow, 'Children');
    
    if ~isempty(plotItems); delete(plotItems); end
        
    set(pH ,...
        'XLim',[0 1],...
        'YLim',[0 1])
    drawnow
    
    % Update GUI workspace
    gui.PTS     = [];
    gui.xyzFun  = [];
    
    % gui.sb.setText( 'Ready')
    
end

% guidata(fig,gui); 

end
