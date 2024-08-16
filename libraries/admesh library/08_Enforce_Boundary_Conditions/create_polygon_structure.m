function POLY = create_polygon_structure(varargin)
% create_polygon_structure - Created polygons around each edge in points
%
% Syntax:  POLY = create_polygon_structure(points,delta)
%
% Inputs:
%    points - points of polygon
%    delta - Grid spacing
%
% Outputs:
%    POLY - Data structure with field Edge 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West, Ethan Kubatko
% The Ohio State University
% email address: dww.425@gmail.com

%------------- BEGIN CODE --------------

%------------------------------------------------------------------------------
% Distribute inputs
%------------------------------------------------------------------------------
% If only one input is given, extend the polygon by the edge length
if size(varargin,2) == 1
    d = [];
    points = varargin{1};
else
    % 
    [points,delta] = deal(varargin{:});
    d = sqrt(3)/2*delta;
end

%------------------------------------------------------------------------------
% Add first 2 coordinates to end of "points" vector to ensure polygon closes
%------------------------------------------------------------------------------
%if isequal(points(1,:), points(end,:)); points(end,:) = []; end
%points = [points; points(1:2,:)];

%------------------------------------------------------------------------------
% Initialize variables
%------------------------------------------------------------------------------
N = size(points,1)-1;
i = 1:N;
dx = zeros(N,1); dy = zeros(N,1);

%------------------------------------------------------------------------------
% Assign P1 P2 P3 (Vectorized)
%------------------------------------------------------------------------------
x1 = points(i,1); x2 = points(i+1,1); 
y1 = points(i,2); y2 = points(i+1,2); 

%------------------------------------------------------------------------------
% Compute the difference between points, edge length, and normal 
%------------------------------------------------------------------------------
dx(:,1) = x1 - x2; 
dy(:,1) = y1 - y2; 
POLY.L = sqrt(dx.^2 + dy.^2);
nx = dy./POLY.L;
ny = -dx./POLY.L;
if isempty(d); d = (sqrt(3)/2).*POLY.L; else d = ones(numel(x1),1).*d;  end

%------------------------------------------------------------------------------
% Create Polygon structure
%------------------------------------------------------------------------------
POLY.x = [ x1 + nx.*d, x1 - nx.*d, x2 - nx.*d,  x2 + nx.*d ];
POLY.y = [ y1 + ny.*d, y1 - ny.*d, y2 - ny.*d,  y2 + ny.*d ];
t = linspace(0,2*pi,500);
%POLY.Circ.x = ones(numel(x1),1)*d*cos(t)+x1*ones(1,500);
%POLY.Circ.y = ones(numel(x1),1)*d*sin(t)+y1*ones(1,500);   
POLY.Circ.x = (ones(numel(x1),1)*cos(t)).*(d*ones(1,500)) + x1*ones(1,500);
POLY.Circ.y = (ones(numel(x1),1)*sin(t)).*(d*ones(1,500)) + y1*ones(1,500);   
% pmidx = (x1+x2)/2;
% pmidy = (y1+y2)/2;  
% POLY.MidCirc.x = ones(numel(pmidx),1)*d*cos(t)+pmidx*ones(1,500);
% POLY.MidCirc.y = ones(numel(pmidy),1)*d*sin(t)+pmidy*ones(1,500);  
end

