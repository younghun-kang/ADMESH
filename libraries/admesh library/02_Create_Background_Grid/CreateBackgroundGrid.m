function [X,Y,delta] = CreateBackgroundGrid(PTS,hmax,hmin,res,UIFigure)
% create_background_mesh - Creates a structured background grid over the
% domain PTS 
%
% Syntax:  [X,Y,delta] = create_background_mesh(PTS,hmax,hmin,res)
%
% Inputs:
%    PTS  - Data structure with fields x & y
%           PTS(1).x = x-coordinates of first polygon
%           PTS(1).y = y-coordinates of first polygon
%    hmin - Minimum element size
%    hmax - Maximum element size
%    res  - a factor multiple, either 1,2,or 3, of hmin controlling the 
%           grid spacing.  
%
% Outputs:
%    X - x-coordinates to rectangular grid
%    Y - y-coordinates to rectangular grid
%    delta - Grid spacing
% Other m-files required: none
% Subfunctions: meshgrid
% MAT-files required: none

% Author: Dustin West, Colton Conroy, Ethan Kubatko
% The Ohio State University
% email address: dww.425@gmail.com
% August 2013; Last revision: 08-August-2013

%------------- BEGIN CODE --------------

%--------------------------------------------------------------------------
% Status bar Update
%--------------------------------------------------------------------------
msg = 'Creating background mesh...';
uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');
%--------------------------------------------------------------------------
% Find boundaing box
%--------------------------------------------------------------------------
xmin = min(vertcat(PTS.Poly(:).x));
xmax = max(vertcat(PTS.Poly(:).x));
ymin = min(vertcat(PTS.Poly(:).y));
ymax = max(vertcat(PTS.Poly(:).y));

%------------------------------------------------------------------------------
% Create cartesian background mesh 
%------------------------------------------------------------------------------
% Determine grid spacing
delta = 1/res*hmin; 

% Create Background Mesh
[X,Y] = meshgrid(xmin-hmax:delta:xmax+hmax,ymin-hmax:delta:ymax+hmax);     

end