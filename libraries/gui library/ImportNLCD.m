function ImportNLCD(app)

if isempty(app.MESH)
    msg = 'No mesh data is found.';
    uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');

    return;
end

MESH = app.MESH;

%--------------------------------------------------------------------------
% Get the file name from the user
%--------------------------------------------------------------------------
app.ProgressBarButton.Text = 'Select a (.tif) file....'; drawnow;

[file, path] = uigetfile({'*.tiff;*.tif','National Land Cover Dataset Files (*.tiff,*.tif)'},...
    'Select a file',cd);

%--------------------------------------------------------------------------
% Did the user make a selection?
%--------------------------------------------------------------------------
if ~file
    app.ProgressBarButton.Text = 'Ready'; drawnow;
    return
end

%--------------------------------------------------------------------------
% Determine the file type
%--------------------------------------------------------------------------
app.ProgressBarButton.Text = 'Reading in data....'; drawnow;

%--------------------------------------------------------------------------
% Read in file
%--------------------------------------------------------------------------
[~,~,ext] = fileparts([path file]);
if any(strcmp(ext,{'.tif','.tiff'}))
    
    filename = [path file];
    Z = imread(filename);
    Z = flipud(Z);
    Z(Z < -1e20) = nan;
    gtinfo = geotiffinfo(filename);

    ax = gtinfo.BoundingBox;
    x = linspace(ax(1,1),ax(2,1),size(Z,2));
    y = linspace(ax(1,2),ax(2,2),size(Z,1));
    [X,Y] = meshgrid(x,y);
    
end

%--------------------------------------------------------------------------
% Retrieve Manning's n value
%--------------------------------------------------------------------------
app.ProgressBarButton.Text = 'Retrieving Manning''s n....'; drawnow;

% Find nearest point
PTS = MESH.Points;
id = knnsearch([X(:),Y(:)],PTS(:,1:2));
NLCD = Z(id);

% Load lookup table
T = LoadNLCDLookupTable('Janssen2016');

% Find indices
[~,id] = ismember(NLCD,T(:,1));
id = nonzeros(id);

% Check if unknown values exist
if length(id) ~= length(NLCD)
    msg = 'Unidentified NLCD values are found.';
    uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');

    return;
end

% Save retrieved Manning's n value
MESH.Attributes.ManningN = T(id,2);

if any([MESH.Constraints.num] == 18 | [MESH.Constraints.num] == 6)
    ix = find([MESH.Constraints.num] == 18 | [MESH.Constraints.num] == 6);
    nodeStr = vertcat(MESH.Constraints(ix).nodeStr);
    nodeStr = unique(nodeStr);

    id = knnsearch([X(:),Y(:)],PTS(nodeStr,1:2));
    NLCD = Z(id);

    % Find indices
    [~,id] = ismember(NLCD,T(:,1));
    id = nonzeros(id);

    MESH.Attributes.ManningN_edge = [double(nodeStr),T(id,2)];
end

app.MESH = MESH;
app.ProgressBarButton.Text = 'Ready'; drawnow;

function T = LoadNLCDLookupTable(ref)
switch lower(ref)
%     case 'bunya2020' % This reference based on 1992 NLCD classification
%         T = [
%             11 0.020 % Open water
%             12 0.022 % Ice/snow
%             21 0.120 % Low residential
%             22 0.121 % High residential
%             23 0.050 % Commercial
%             31 0.040 % Bare rock/sand
%             32 0.060 % Gravel pit
%             33 0.100 % Transitional
%             41 0.160 % Deciduous forest
%             42 0.180 % Evergreen forest
%             43 0.170 % Mixed forest
%             51 0.070 % Shrub land
%             61 0.100 % Orchard/vineyard
%             71 0.035 % Grassland
%             81 0.033 % Pasture
%             82 0.040 % Row crops
%             83 0.035 % Small grains
%             84 0.032 % Fallow
%             85 0.030 % Recreational grass
%             91 0.140 % Woody wetland
%             92 0.035 % Herbaceous wetland
%             ];
    case 'kalyanapu2009'
        T = [
            21 0.0404 % Developed, open space
            22 0.0678 % Developed, low intensity
            23 0.0678 % Developed, medium intensity
            24 0.0404 % Developed, high intensity
            31 0.0113 % Barren land
            41 0.36   % Deciduous forest
            42 0.32   % Evergreen forest
            43 0.40   % Mixed forest
            52 0.40   % Shrub/scrub
            71 0.368  % Grassland/herbaceous
            81 0.325  % Pasture/Hay
            90 0.086  % Woody wetlands
            95 0.1825 % Emergent herbaceous wetlands
            ];
    case 'janssen2016'
        T = [
            11 0.040 % Open water
            21 0.040 % Developed, open space 
            22 0.100 % Developed, low intensity
            23 0.080 % Developed, medium intensity
            24 0.150 % Developed, high intensity
            31 0.025 % Barren land
            41 0.160 % Deciduous forest
            42 0.160 % Evergreen forest
            43 0.160 % Mixed forest
            52 0.100 % Shrub/scrub
            71 0.035 % Grassland/herbaceous
            81 0.030 % Pasture/Hay
            82 0.035 % Cultivated crops
            90 0.120 % Woody wetlands
            95 0.070 % Emergent herbaceous wetlands
            ];
end

end


end







