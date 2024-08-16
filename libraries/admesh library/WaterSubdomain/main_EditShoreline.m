addpath(genpath('admesh-lib-part'));
%--------------------------------------------------------------------------
% Load example file (variable "Points")
%--------------------------------------------------------------------------
ShorelineFile = ['D:\academic\data\pr-preevents\lower_Neches\ArcGIS\',...
                'ClintExtent/Shoreline_InBoundary.shp'];
Shoreline     = shaperead(ShorelineFile);

for i = 1 : length(Shoreline)
    Shoreline(i).np = length(Shoreline(i).X);
end

close all;

fid = 0 : length(Shoreline)-1;

id = fid+1;

pgon_x = []; pgon_y = [];
for i = 1 : length(id)
    j = id(i);
    pgon_x{i} = Shoreline(j).X;
    pgon_y{i} = Shoreline(j).Y;
end
pgon = polyshape(pgon_x,pgon_y);
pgon = rmholes(pgon);
%%
% plot(pgon_edit); axis equal;
roi = drawpolygon; roi = polyshape(roi.Position);
pgon_edit = subtract(pgon_edit,roi);


I = find(isnan(pgon_edit.Vertices(:,1)));
I(:,2) = I(:,1)+1;
I = I';
I = [1; I(:); size(pgon_edit.Vertices,1)+1];
I = diff(I);
temp = mat2cell(pgon_edit.Vertices,I,2);
temp = temp(1:2:end);

I = cellfun(@(x) length(x) < 50,temp);
temp(I) = [];

pgon_x = cellfun(@(x) x(:,1),temp,'UniformOutput',0);
pgon_y = cellfun(@(x) x(:,2),temp,'UniformOutput',0);
pgon_edit2 = polyshape(pgon_x,pgon_y);
hold on; plot(pgon_edit2);




