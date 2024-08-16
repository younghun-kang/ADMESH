function SaveMatFile(varargin)
% SaveMatFile - GUI Callback that saves PTS data structure 
%
% Syntax:  SaveMatFile(varargin)
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

%--------------------------- BEGIN CODE -----------------------------------

%--------------------------------------------------------------------------
% Get GUI data
%--------------------------------------------------------------------------
app = varargin{1};

if isempty(app.PTS) % User has not run ADmesh yet
    msg = 'No edge structure to save....';
    uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');
    return
end

% Save edge structure
% gui.sb.setText('Saving data in a (.mat) file...')

[filename,pathname] = uiputfile('*.mat','Save File As');

if filename ~=0
    
    PTS = app.PTS; %#ok<*NASGU>
    xyzFun = app.xyzFun;
    
    % Save data to mat file
    save([pathname, filename], 'PTS', 'xyzFun', '-v7.3');
    
    app.FilePath = [pathname filename]; % Save current file path
    
    % guidata(fig,gui); % Update GUI data
    
    % gui.sb.setText('Ready')
    
else
    % gui.sb.setText('Ready')
    return
end