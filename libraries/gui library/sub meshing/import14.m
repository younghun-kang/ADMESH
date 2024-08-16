function [PTS,status,xyzIncluded,xyzFun] = import14(varargin)

xyzIncluded = 0;
status = 1;
PTS = [];

%-----------------------------------------------------------------------
% Get filename & location from user
%-----------------------------------------------------------------------
uiStatusBar('Select a file...')

%   Get the file name from the user
[filename, pathname] = uigetfile(...
    {'*.14;*.grd;*.2dm','Grid Files (*.14,*.grd,*.2dm)'},'Select a file');

% If user cancels
if filename == 0
    
    uiStatusBar('Ready')
    
    status = 0;
    
    return
    
end

%------------------------------------------------------------------------------
% Assign filename & path to guiH structure & read in file based on file
% extension
%------------------------------------------------------------------------------

% Assign filename & path to guiH data structure
FilePath = [pathname filename]; clear pathname filename

% Determine the file type
[~,~,ext] = fileparts(FilePath);

if( strcmp(ext,'.14') || strcmp(ext,'.grd') ) % Read in ADCIRC file
    
    mesh = Read14File(FilePath);
    
elseif strcmp(ext,'.2dm') % Read in general SMS file
    
    mesh = READ_2DM(FilePath);
    
end

if isempty(mesh); status = 0; return; end

[PTS,xyzFun] = Mesh2PTS(mesh);

if isa(xyzFun,'griddedInterpolant')
    
    msg = 'What coordinate system would you like to write your data in?';
    choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'Yes','No'},'DefaultOption',2,'Icon','Warning');
    switch choice
        
        case 'Yes'
            xyzIncluded = 1;
        case 'No'
            xyzFun = [];
    end
    
end













