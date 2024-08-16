function [PTS,status,xyzIncluded,xyzFun] = importKMLv2(varargin)

xyzIncluded = 0;
xyzFun = [];
status = 1;
PTS = [];

%---------------------------------------------------------------------
% Read in KML file
%---------------------------------------------------------------------
% Ask user for file info
[filename,pathname] = uigetfile('*.kml','Select KML file...');

% If user cancels
if filename == 0
    
    uiStatusBar('Ready')
    
    status = 0;
    
    return
    
end

% Open file
fid = fopen([pathname filename],'r');

uiStatusBar('Reading in KML file...')

% Read file into cell
C = textscan(fid,'%s');

% Close file id
fclose(fid);

%---------------------------------------------------------------------
% What kind of geometry are we reading in?
%---------------------------------------------------------------------
polygon = ~isempty(find(~cellfun(@isempty,strfind(C{1}, '<Polygon>')), 1));

% Check if file contains polygons and/or lines
if polygon == 0
    
    msg = 'To be imported a KML file should contain a polygon.';
    uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');
    uiStatusBar('Ready')
    status = 0;
    return
    
end

% Search for beginning coordinate flags
idx = find(~cellfun(@isempty,strfind(C{1}, '<coordinates>')))+1;

% Search for ending coordinate flags
idy = find(~cellfun(@isempty,strfind(C{1}, '</coordinates>')))-1;

% Initialize cell for holding segments
S = cell(numel(idx),1);

% Initialize area vector
A = zeros(numel(S),1);

% Store path segments in cell

uiStatusBar('Compiling polygon(s)...')

for k = 1:numel(idx)
    
    % Extract section
    xyz =  cell2mat(cellfun(@str2num,C{1}(idx(k):idy(k)),'uniformoutput',0));
    
    % Convert to cartesian coordinates
    if k == 1
        [x,y,PTS.cpplon,PTS.cpplat] = Geo2Meters(xyz(:,1),xyz(:,2));
    else
        [x,y] = Geo2Meters(xyz(:,1),xyz(:,2),PTS.cpplon,PTS.cpplat);
    end
    
    % Make sure point set is unique and closed
    S{k} = [unique([x y],'rows','stable'); [x(1) y(1)]];
    
    % Compute polygon area
    A(k) = polyarea(S{k}(:,1),S{k}(:,2));
    
    uiStatusBar(k/numel(idx))
    
end

% Sorting polygons based on area
[~,i] = sort(A,'descend');

% Create edge structure

for k = 1:numel(S)
    
    PTS.Poly(k).x = S{i(k)}(:,1);
    PTS.Poly(k).y = S{i(k)}(:,2);
    
    uiStatusBar(k/numel(S))
    
end



end