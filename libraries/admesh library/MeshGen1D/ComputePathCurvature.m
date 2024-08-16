function [PI,p] = ComputePathCurvature(XY,FixedPointsXY,SmoothingRMSE,varargin)
%==========================================================================
% ComputeCurveCurvature_v1
% - This function computes curvature and total length of given curve.
% - Outputs 
%   PI        : Projection of xy plane to unit inverval [0,1].
%               This includes projections of positions, first and second order derivatives.
%   CurvLength: Length of input curve which is numerically computed
%   p         : Coordinates in 1D domain [0,1] which corresponds to input xy
% 
% Update history
% (v1) 2021-04-07 Written by Younghun Kang
% (v1) 2021-04-11 YK: Pass projection functions as one output
% 
%==========================================================================

%==========================================================================
% Parser optional arguments
%==========================================================================
%--------------------------------------------------------------------------
% Set default values
%--------------------------------------------------------------------------
WARN = 'on';
%--------------------------------------------------------------------------
% Get optional values
%--------------------------------------------------------------------------
SkipNext = false;
for i = 1 : length(varargin)
    if SkipNext
        SkipNext = false;
        continue;
    end
    ivar = varargin{i};
    if strcmpi(ivar,'warning')
        WARN = varargin{i+1};
        SkipNext = true;
    else
        error('Incorrect use of optional arguments.');
    end
end

%--------------------------------------------------------------------------
% Find index of fixed points
%--------------------------------------------------------------------------
id_fixedPoints = [];
if ~isempty(FixedPointsXY)
    id_fixedPoints = find(ismember(XY,FixedPointsXY,'rows'));
end

%--------------------------------------------------------------------------
% Apply smoothing
%--------------------------------------------------------------------------
x0 = XY(:,1); y0 = XY(:,2);
if isequal(XY(1,:),XY(end,:)) % if periodic
    x0 = x0(1:end-1);
    y0 = y0(1:end-1);
    x0 = repmat(x0,3,1);
    y0 = repmat(y0,3,1);
    
    id_fixedPoints(id_fixedPoints == size(XY,1)) = [];
    n = size(XY,1) - 1;
    id_fixedPoints = [id_fixedPoints; n+id_fixedPoints; 2*n+id_fixedPoints];
end
%---1D coordinate vector of raw path (not sure how it should be used)
p1 = cumsum(sqrt(sum(diff(XY).^2,2)));
p1 = [0; p1(:)];
if length(p1) < 3
    PI = [];
    return;
end
if SmoothingRMSE >= 0
%     t = 1 : length(x);
%     Smooth_x = fit(t(:),x(:),'smoothingspline','SmoothingParam',Setting.SmoothingParam);
%     Smooth_y = fit(t(:),y(:),'smoothingspline','SmoothingParam',Setting.SmoothingParam);
%     x = Smooth_x(t);
%     y = Smooth_y(t);
    
    %---enforce fixed points
    t = 1 : length(x0);
    w = ones(length(x0),1);
    w(id_fixedPoints) = 1e8;
    
    SmoothingParam = 0;
    Smooth_x = csaps(t,x0,SmoothingParam,t,w);
    Smooth_x(id_fixedPoints) = x0(id_fixedPoints);
    
    Smooth_y = csaps(t,y0,SmoothingParam,t,w);
    Smooth_y(id_fixedPoints) = y0(id_fixedPoints);
    RMSE = sqrt( sum( (x0(:) - Smooth_x(:)).^2 + (y0(:) - Smooth_y(:)).^2)/length(x0));
    if RMSE < SmoothingRMSE
        SmoothingParam = 1e-10;
        if strcmpi(WARN,'on')
        warning(['The input RMSE cannot be reached for this segment. ',...
            'Skip iteration and result straight line segment (smoothing param = 0).']);
        end
    else
    SP1 = 0;
    SP2 = 1;
    while 1
        if SP2 - SP1 < 1e-16
            if strcmpi(WARN,'on')
            warning(['Iteration cannot reach to the stopping criterion. ',...
                'It is stopped and result with RMSE closest to the given value.']);
            end
            break;
        end
        SmoothingParam = (SP1 + SP2)/2;
        
        Smooth_x = csaps(t,x0,SmoothingParam,t,w);
        Smooth_x(id_fixedPoints) = x0(id_fixedPoints);
        
        Smooth_y = csaps(t,y0,SmoothingParam,t,w);
        Smooth_y(id_fixedPoints) = y0(id_fixedPoints);
        RMSE = sqrt( sum( (x0(:) - Smooth_x(:)).^2 + (y0(:) - Smooth_y(:)).^2)/length(x0));
        if RMSE < SmoothingRMSE
            SP2 = SmoothingParam;
        else
            SP1 = SmoothingParam;
        end
        
        if abs(RMSE - SmoothingRMSE)/SmoothingRMSE < 1e-3
            break;
        end
    end
    end
    x = Smooth_x;
    y = Smooth_y;
else
    x = x0;
    y = y0;
end

%--------------------------------------------------------------------------
% Compute spline and derivative of spline
%--------------------------------------------------------------------------
t = 1 : length(x);
sx_t = csape(t,x);
sy_t = csape(t,y);
dsx_t  = @(t) ppval(fnder(sx_t,1),t);
dsy_t  = @(t) ppval(fnder(sy_t,1),t);

%--------------------------------------------------------------------------
% Compute distance
%--------------------------------------------------------------------------
% farc_t = @(t) sqrt( dsx_t(t).^2 + dsy_t(t).^2);
% CurvDist = zeros(length(t),1);
% CurvDist(1) = 0;
% for i = 2 : length(t)
%     CurvDist(i) = CurvDist(i-1) + integral(farc_t,t(i-1),t(i),'reltol',1e-10);
% end

%--------------------------------------------------------------------------
% Compute distance using trapezoidal rule
%--------------------------------------------------------------------------
farc_t = @(t) sqrt( dsx_t(t).^2 + dsy_t(t).^2);
t1 = t(1:end-1);
t2 = t(2:end);
t3 = (t1+t2)/2;
% CurvDist = (farc_t(t1) + farc_t(t2))/2;
CurvDist = (farc_t(t1) + farc_t(t2) + 2*farc_t(t3))/4;
CurvDist = cumsum(CurvDist);
CurvDist = [0; CurvDist(:)];

%--------------------------------------------------------------------------
% Normalize distance vector so that 1D domain becomes unit interval (0,1)
%--------------------------------------------------------------------------
p = CurvDist;
% p = p*t(end)/p(end);
if isequal(XY(1,:),XY(end,:)) % if periodic
    if isempty(id_fixedPoints)
        m = n+1;
    else
        m = n+id_fixedPoints(1);
    end
    p = p - p(m);
end
%--------------------------------------------------------------------------
% Cubic spline within 1D coordinate system
%--------------------------------------------------------------------------
sx   = csaps(t,x0,SmoothingParam,[],w);
sy   = csaps(t,y0,SmoothingParam,[],w);
p2t = @(s) interp1(p,t,s);
dsx  = @(s) ppval(fnder(sx,1),p2t(s));
ddsx = @(s) ppval(fnder(sx,2),p2t(s));
dsy  = @(s) ppval(fnder(sy,1),p2t(s));
ddsy = @(s) ppval(fnder(sy,2),p2t(s));

%--------------------------------------------------------------------------
% Compute curvature
%--------------------------------------------------------------------------
kappa = @(s) abs((dsx(s).*ddsy(s) - dsy(s).*ddsx(s))./(dsx(s).^2 + dsy(s).^2).^(3/2));

if any(isnan(kappa(p)))
    0;
end
%--------------------------------------------------------------------------
% Store projection functions for output
%--------------------------------------------------------------------------
PI.x = @(s) interp1(p,x0,s);
PI.y = @(s) interp1(p,y0,s);
PI.x1 = PI.x;
PI.y1 = PI.y;
PI.sx = @(s) interp1(p,x,s);
PI.sy = @(s) interp1(p,y,s);
PI.k  = kappa;

if isequal(XY(1,:),XY(end,:)) % if periodic
    p = p(m+(0:n));
end
PI.p  = p;

PI.dx  = dsx ;
PI.ddx = ddsx;
PI.dy  = dsy ;
PI.ddy = ddsy;

PI.SmoothingParam = SmoothingParam;

