function [PTS,status,xyzIncluded,xyzFun] = importMAT(varargin)

xyzIncluded = 0;
status = 1;
PTS = [];

%-----------------------------------------------------------------------
% Get filename & location from user
%-----------------------------------------------------------------------
uiStatusBar('Select a file...')

%   Get the file name from the user
[filename, pathname] = uigetfile({'*.mat'},'Select a file');

% If user cancels
if filename == 0
    
    uiStatusBar('Ready')
    
    status = 0;
    
    return
    
end

%-----------------------------------------------------------------------
% Assign filename & path to guiH structure
%-----------------------------------------------------------------------
FilePath = [pathname filename];

%-----------------------------------------------------------------------
% Check file
%-----------------------------------------------------------------------
uiStatusBar('Checking file...')

matData = struct2cell(whos('-file',[pathname filename]));

if ~any(strcmp(matData(1,:),'PTS'))
    
    msg = ['No edge structure exists in this file.' ...
        ' Make sure you are selecting the correct file for ADMESH.'];
    uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');
    
    uiStatusBar('Ready')
    
    status = 0;
    
    return
    
end

%-----------------------------------------------------------------------
% Load edge structure
%-----------------------------------------------------------------------
uiStatusBar('Loading edge structure...')


if any(strcmp(matData(1,:),'PTS'))
    
    % Initialize as 0
    PTS = 0;
    
    % Load PTS variable
    load(FilePath, 'PTS')
    
    % Re-format to new ADMESH standard
    if ~isfield(PTS,'Poly')
        
        temp = PTS; clear PTS
        
        PTS.Poly(numel(temp),1) = struct('x',[],'y',[]);
        
        for i = 1:numel(temp)
            PTS.Poly(i).x = temp(i).x;
            PTS.Poly(i).y = temp(i).y;
        end
        
        clear temp
        
    end
    
end

%-----------------------------------------------------------------------
% Load elevation data
%-----------------------------------------------------------------------
uiStatusBar('Loading elevation data...')

if any(strcmp(matData(1,:),'xyzFun'))
    
    % Initialize as 0
    xyzFun = [];
    
    % Load xyzFun variable
    load(FilePath, 'xyzFun')
    
    if isa(xyzFun,'griddedInterpolant')
        
        msg = ['Elevation data was detected in this file. Would you ',...
            'prefer to use this data set?'];
        choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
            'Options',{'Yes','No'},'DefaultOption',2,'Icon','Warning');
        switch choice
            
            case 'Yes'
                xyzIncluded = 1;
            case 'No'
                xyzFun = [];
        end
        
    end
    
end

%------------------------------------------------------------------------------
% Load geo coordinate references. If they exist convert to new
%------------------------------------------------------------------------------
if any(strcmp(matData(1,:),'CPPLAT'))
    
    load(guiH.FilePath, 'CPPLAT')
    
    PTS.cpplat = CPPLAT;
    
end

if any(strcmp(matData(1,:),'CPPLON'))
    
    load(guiH.FilePath, 'CPPLON')
    
    PTS.cpplon = CPPLON;
    
end


%------------------------------------------------------------------------------
% Check for NaN's in cpplon and cpplat and remove in accordance with new
% Admesh standard
%------------------------------------------------------------------------------
if isfield(PTS,'cpplon')
    if isnan(PTS.cpplon)
        PTS = rmfield(PTS,'cpplon');
        PTS = rmfield(PTS,'cpplat');
    end
end

end