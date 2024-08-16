function ImportMesh(varargin)
% OpenMeshCallback - GUI Callback that calls on Read_14 to read in a fort.14
% file.
%
% Syntax:  Open14Callback(guiFig,~,~)
%
% Inputs:
%    guiFig - handle that identifies the figure
%
% Outputs:
%    non
%
% Other m-files required: READ_14
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
%---------------------------- BEGIN CODE --------------------------------------

%------------------------------------------------------------------------------
% Get GUI data
%------------------------------------------------------------------------------
if nargin == 1
    
    app = varargin{1};
    
    msg = 'Select a file...';
    progdlg = uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');
    
    % Ask the user to select a file
    [filename, pathname] = uigetfile(...
        {'*.14;*.grd;*.2dm','Grid Files (*.14,*.grd,*.2dm)'},'Select a file');
    
    if filename == 0
        return
    end
    
    % Assign filename & path to guiH data structure
    app.FilePath = [pathname filename]; clear pathname filename
    
else
        
end

% Turn off colormap
SetContourStatus(app,'Off');

%------------------------------------------------------------------------------
% Erase any sub meshing stuff
%------------------------------------------------------------------------------
CancelADMESHSubMesh(fig,nan,'c')

%-------------------------------------------------------------------------
% Clear guiH data structure variables
%-------------------------------------------------------------------------
app.PTS        = [];
app.xyzFun     = [];
app.MESH       = [];

%------------------------------------------------------------------------------
% Assign filename & path to guiH structure & read in file based on file
% extension
%------------------------------------------------------------------------------
% Determine the file type
[~,~,ext] = fileparts(app.FilePath);

if( strcmp(ext,'.14') || strcmp(ext,'.grd') ) % Read in ADCIRC file
    
    [MESH,xyzFun,status] = Read14File(app.FilePath,app);
    
    if status == 0
        % gui.sb.setText('Ready.');
        errordlg('There was an error reading in the file.','ADMESH')
        return
    end
    
elseif strcmp(ext,'.2dm') % Read in general SMS file
    
    %mesh = READ_2DM(gui.FilePath);
    
end


% gui.sb.setText('Creating ADMESH edge structure...');
if ~isempty(MESH.Constraints)
    
    % Check for external boundary constraints
    I = find(ismember([MESH.Constraints.num],[0 2 10 12 20 22 30 -1 3 13 23]), 1);
    
    if ~isempty(I)
        PTS = Constraints2PTS(MESH);
    else
        PTS = triangulation2PTS(MESH);
    end
else
    
    PTS = triangulation2PTS(MESH);

end

% Save edge structure
% gui.sb.setText('Saving data in a (.mat) file...');

[filename,pathname] = uiputfile('*.mat','Save File As');

if filename ~=0
    
    % Save data to mat file
    save([pathname, filename], 'PTS', 'xyzFun', '-v7.3');
    
    app.FilePath = [pathname filename]; % Save current file path
    app.PTS      = PTS;
    app.xyzFun   = xyzFun;
    
    % guidata(fig,gui); % Update GUI data
    
    clear guiH mesh
    
    % Plot Boundary
    PlotEdgeStructure(app,.1);
    
    drawnow
    
else
    % gui.sb.setText('Ready');
    return
end

% gui.sb.setText('Ready');

end