function h0 = BathymetryFunction(h0,X,Y,Z,s,hmin,hmax,delta,Settings,UIFigure)
% bathymetry - Computes the mesh size based on bathymetry changes
%
% Syntax:  [h_bathy, B_grid] = bathymetry(X,Y,D,bathy,bathy_value,delta)
%
% Inputs:
%    X - x-coordinates to rectangular grid
%    Y - y-coordinates to rectangular grid
%    D - Distance function
%    xyzFun - bathymetry data
%    s - bathymetry mesh size factor
%    delta - background grid spacing
%    guiFig - handle that identifies the figure
%
% Outputs:
%    h_bathy - mesh size based on bathymetry
%    Z - bathymetry interpolated onto background grid
%
% Other m-files required: none
% Subfunctions: InPolygon,bwdist,Compute_Distance_v3
% MAT-files required: none
%
% Author: Dustin West, Colton Conroy, Ethan Kubatko
% The Ohio State University
% email address: dww.425@gmail.com,conroy.51@osu.edu
% August 2013; Last revision: 08-August-2013

%------------- BEGIN CODE --------------

%------------------------------------------------------------------------------
% Elevation Status
%------------------------------------------------------------------------------
Status = Settings.B.Status;

if strcmp(Status,'On') && ~isempty(Z) % If true, compute bathymetry function
    
    %------------------------------------------------------------------------------
    % Initialize variables & status
    %------------------------------------------------------------------------------
    twodeltax = 1/(2*delta);
    twodeltay = 1/(2*delta);
    gradBx = zeros(size(X));
    gradBy = zeros(size(Y));
    [LY,LX] = size(X);
        
    %------------------------------------------------------------------------------
    % Compute the gradient of the bathymetry
    %------------------------------------------------------------------------------
    msg = 'Computing the gradient of the bathymetry...';
    uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');

    % Inner boundary
    i = (2:LX-1); j = (2:LY-1);
    gradBx(j,i) = (Z(j,i+1)-Z(j,i-1))*twodeltax;
    gradBy(j,i) = (Z(j+1,i)-Z(j-1,i))*twodeltay;
    
    % Left boundary
    i = 1; j = (2:LY-1);
    gradBx(j,i) = (-3*Z(j,i)+4*Z(j,i+1)-Z(j,i+2))*twodeltax;
    gradBy(j,i) = (Z(j+1,i)-Z(j-1,i))*twodeltay;
    
    % Right boundary
    i = LX; j = 2:LY-1;
    gradBx(j,i) = (Z(j,i-2)-4*Z(j,i-1)+3*Z(j,i))*twodeltax;
    gradBy(j,i) = (Z(j+1,i)-Z(j-1,i))*twodeltay;
    
    
    % Top boundary
    j = LY; i = (2:LX-1);
    gradBx(j,i) = (Z(j,i+1)-Z(j,i-1))*twodeltax;
    gradBy(j,i) = (Z(j-2,i)-4*Z(j-1,i)+3*Z(j,i))*twodeltay;
    
    
    % Bottom boundary
    j = 1; i = (2:LX-1);
    gradBx(j,i) = (Z(j,i+1)-Z(j,i-1))*twodeltax;
    gradBy(j,i) = (-3*Z(j,i)+4*Z(j+1,i)-Z(j+2,i))*twodeltay;
    
    
    % Top left corner
    i = 1;
    j = LY;
    gradBx(j,i) = (-3*Z(j,i)+4*Z(j,i+1)-Z(j,i+2))*twodeltax;
    gradBy(j,i) = (Z(j-2,i)-4*Z(j-1,i)+3*Z(j,i))*twodeltay;
    
    % Bottom left corner
    i = 1;
    j = 1;
    gradBx(j,i) = (-3*Z(j,i)+4*Z(i+1,j)-Z(i+2,j))*twodeltax;
    gradBy(j,i) = (-3*Z(j,i)+4*Z(j+1,i)-Z(j+2,i))*twodeltay;
    
    % Bottom right corner
    i = LX;
    j = 1;
    gradBx(j,i) = (Z(j,i-2)-4*Z(j,i-1)+3*Z(j,i))*twodeltax;
    gradBy(j,i) = (-3*Z(j,i)+4*Z(j+1,i)-Z(j+2,i))*twodeltay;
    
    % Top right corner
    i = LX;
    j = LY;
    gradBx(j,i) = (Z(j,i-2)-4*Z(j,i-1)+3*Z(j,i))*twodeltax;
    gradBy(j,i) = (Z(j-2,i)-4*Z(j-1,i)+3*Z(j,i))*twodeltay;
    
    msg = 'Computing mesh size from bathymetry...';
    uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');

    % Calculate hb(x), the bathymetry mesh size
    h_bathy = s.*(abs(sqrt((gradBx).^2 + (gradBy).^2)./Z)).^(-1);
    
%     if min(h_bathy(:)) < hmin % Normalize
% 
%         % Find absolute min & max values
%         %minVal = min(h_bathy(:));
%         %maxVal = max(h_bathy(:));
%         minVal = hmin;
%         maxVal = hmax;
%         
%         % Normalize
%         h_bathy = abs(h_bathy - minVal) ./ ( maxVal - minVal );
%         
%     end
    
    % A
    
    h_bathy(h_bathy == 0) = hmax;
    
    % Enforce boundary conditions
    h_bathy(h_bathy < hmin) = hmin;
    h_bathy(h_bathy > hmax) = hmax;
    
    %figure
    %contourf(X,Y,h_bathy)
    %colorbar
 
    % Compare initial conditions and save
    h0 = min(h_bathy, h0);
        
end
end