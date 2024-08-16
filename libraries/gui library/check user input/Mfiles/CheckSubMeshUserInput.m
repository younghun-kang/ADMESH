function quit = CheckSubMeshUserInput(fig)
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
gui = % guidata(fig);

% Check minimum element size
if isempty(strtrim(get(gui.MinElementSize,'string')));
    msg = 'Enter a minimum element size'; 
    uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');
    quit = 1;
    return
end

% Check maximum element size
if isempty(strtrim(get(gui.MaxElementSize,'string')));
    msg = 'Enter a maximum element size';
    uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');
    quit = 1;
    return
end

% Check curvature value if curvature is on
if get(gui.CurvatureStatus, 'value') == 1
    if isempty(strtrim(get(gui.CurvatureValue,'string')));
        msg = 'Enter a curvature value';
        uiconfirm(app.UIFigure,msg,'ADMESH',...
            'Options',{'OK'},'DefaultOption',1,'Icon','Error');
        quit = 1;
        return
    end
end

% Check lfs value if lfs is on
if get(gui.LFSStatus,'value') == 2
    if isempty(strtrim(get(gui.LFSValue ,'string')));
        msg = 'Enter a local feature size value';
        uiconfirm(app.UIFigure,msg,'ADMESH',...
            'Options',{'OK'},'DefaultOption',1,'Icon','Error');
        quit = 1;
        return
    end
end

% Check bathymetry data if bathymetry is on
if get(gui.ElevStatus,'value') == 2
    if isempty(strtrim(get(gui.ElevValue,'string')));
        msg = 'Enter a bathymetry/topography value';
        uiconfirm(app.UIFigure,msg,'ADMESH',...
            'Options',{'OK'},'DefaultOption',1,'Icon','Error');
        quit = 1;
        return
    end
    if ~isa(gui.xyzFun, 'griddedInterpolant');
        msg = ['You must load bathymetry data '...
            'if you would like to use the bathymetry parameter'];
        uiconfirm(app.UIFigure,msg,'ADMESH',...
            'Options',{'OK'},'DefaultOption',1,'Icon','Error');
        quit = 1;
        return
    end
end

% Check dominate tide value if dominate tide is on
if get(gui.TidalStatus,'value') > 1
    if isempty(strtrim(get(gui.TidalValue,'string')));
        msg = 'Enter a dominate tide value';
        uiconfirm(app.UIFigure,msg,'ADMESH',...
            'Options',{'OK'},'DefaultOption',1,'Icon','Error');
        quit = 1;
        return
    end
    if ~(isa(gui.xyzFun, 'griddedInterpolant'));
        msg = ['You must load bathymetry data '...
            'if you would like to use the bathymetry parameter'];
        uiconfirm(app.UIFigure,msg,'ADMESH',...
            'Options',{'OK'},'DefaultOption',1,'Icon','Error');
        quit = 1;
        return
    end
end

% Check mesh grading value if mesh grading is on
if get(gui.GradingStatus,'value') == 1
    if isempty(strtrim(get(gui.GradingValue, 'string')));
        msg = 'Enter a mesh grading value';
        uiconfirm(app.UIFigure,msg,'ADMESH',...
            'Options',{'OK'},'DefaultOption',1,'Icon','Error');
        quit = 1;
        return
    end
end


end