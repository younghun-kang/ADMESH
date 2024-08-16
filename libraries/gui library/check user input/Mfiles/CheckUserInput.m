function quit = CheckUserInput(app)
% CheckUserInput - Checks user input
%
% Syntax:  quit = CheckUserInput(varargin)
%
% Inputs:
%    None
%
% Outputs:
%    quit - 1 to stop the program. 0 means everything is cool
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% Last revision: 26-April-2014

%---------------------------------------------------------------------
% Begin Code
%---------------------------------------------------------------------

% Initialize output
quit = 0;

% Get GUI data
% gui = guidata(fig);

% 1) Check if any data isloaded
if isempty(app.PTS) && isempty(app.MESH)
    
    errordlg('No file is currently loaded.','Error');
    
    quit = 1;
    
    return
    
end

% 2) Check for mesh and no edge structure
if isempty(app.PTS) && ~isempty(app.MESH)
    
    msgbox([...
        'ADMESH needs to import your current mesh into an edge structure '...
        'file. Import your mesh file by selecting Import ---> Mesh.'],'ADMESH')
    
    drawnow; pause(.005);
    
    return    
end

% 3) Check for mesh in plot window
% h = findobj(gui.ViewAxes,'Tag','Mesh'); % Mesh plot tag

% Check for mesh in plot window
if ~isempty(app.MESH)

    % Ask user if they are sure they want to delete the mesh
    msg = 'Are you sure you want to clear this mesh?';
    choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'Clear','Cancel'},'DefaultOption',2,'Icon','Warning');
    
    drawnow; pause(0.05);  % this innocent line prevents the Matlab hang
    
    if strcmp(choice,'Cancel')
        % gui.sb.setText('Ready');
        quit = 1;
        return
    end

    % Clear guiH mesh data structure variables
    app.MESH  = [];

    % Clear patch data
%     delete(h)
    colorbar('delete') % Delete current colorbar

    % Update GUI data
    % guidata(fig,gui); % Update GUI data
    
    % Clear mesh from screen & Plot Boundary
    PlotEdgeStructure(app,.1);
    
end


% 4) Check ADMESH parameters

% Check minimum element size
if app.MinElementSizeEditField.Value <= 0
    msg = 'Enter a minimum element size';
    uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');
    quit = 1;
    return
end

% Check maximum element size
if app.MaxElementSizeEditField.Value <= 0
    msg = 'Enter a maximum element size';
    uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');
    quit = 1;
    return
end

% Check curvature value if curvature is on
if strcmpi(app.BoundaryCurvatureDropDown.Value, 'on')
    if app.BoundaryCurvatureEditField.Value <= 0
        msg = 'Enter a curvature value';
        uiconfirm(app.UIFigure,msg,'ADMESH',...
            'Options',{'OK'},'DefaultOption',1,'Icon','Error');
        quit = 1;
        return
    end
end

% Check lfs value if lfs is on
if strcmpi(app.LocalFeatureSizeDropDown.Value, 'on')
    if app.LocalFeatureSizeEditField.Value <= 0
        msg = 'Enter a local feature size value';
        uiconfirm(app.UIFigure,msg,'ADMESH',...
            'Options',{'OK'},'DefaultOption',1,'Icon','Error');
        quit = 1;
        return
    end
end

% Check bathymetry data if bathymetry is on
if strcmpi(app.ElevationGradientsDropDown.Value, 'on')
    if app.ElevationGradientsEditField.Value <= 0
        msg = 'Enter a bathymetry/topography value';
        uiconfirm(app.UIFigure,msg,'ADMESH',...
            'Options',{'OK'},'DefaultOption',1,'Icon','Error');
        quit = 1;
        return
    end
    if isempty(app.xyzFun)
        msg = ['You must load bathymetry data '...
            'if you would like to use the bathymetry parameter'];
        uiconfirm(app.UIFigure,msg,'ADMESH',...
            'Options',{'OK'},'DefaultOption',1,'Icon','Error');
        quit = 1;
        return
    end
end

% Check dominate tide value if dominate tide is on
if strcmpi(app.DominateTideDropDown.Value, 'on')
    if app.DominateTideEditField.Value <= 0
        msg = 'Enter a dominate tide value';
        uiconfirm(app.UIFigure,msg,'ADMESH',...
            'Options',{'OK'},'DefaultOption',1,'Icon','Error');
        quit = 1;
        return
    end
    if isempty(app.xyzFun)
        msg = ['You must load bathymetry data '...
            'if you would like to use the bathymetry parameter'];
        uiconfirm(app.UIFigure,msg,'ADMESH',...
            'Options',{'OK'},'DefaultOption',1,'Icon','Error');
        quit = 1;
        return
    end
end

% Check mesh grading value if mesh grading is on
if strcmpi(app.MeshGradingDropDown.Value, 'on')
    if app.MeshGradingEditField.Value <= 0
        msg = 'Enter a mesh grading value';
        uiconfirm(app.UIFigure,msg,'ADMESH',...
            'Options',{'OK'},'DefaultOption',1,'Icon','Error');
        quit = 1;
        return
    end
end

end