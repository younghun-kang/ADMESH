function Mesh1D = MeshGeneration1D(PI,fixedPoints,Settings,h2d)
%--------------------------------------------------------------------------
% Get input parameters
%--------------------------------------------------------------------------
h_min = Settings.hmin;
h_max = Settings.hmax;
K     = Settings.K.Value;
g     = Settings.G.Value;

%--------------------------------------------------------------------------
% Loop over 1D projection map array
%--------------------------------------------------------------------------
x = PI.x(PI.p);
y = PI.y(PI.p);

id = ismember([x,y],fixedPoints,'rows');
fixedPoints1D = PI.p(id);
    
%--------------------------------------------------------------------------
% Exit the program if the length of flowline is shorter than h_max
%--------------------------------------------------------------------------
if PI.p(end) < h_min
    Mesh1D = [];
    Mesh1D.X = PI.x(PI.p([1 end]));
    Mesh1D.Y = PI.y(PI.p([1 end]));
    Mesh1D.K = [0; 0];
    Mesh1D.p = PI.p([1 end]);
    Mesh1D.h = [];
    
    return;
end

%--------------------------------------------------------------------------
% Set desired element size
%--------------------------------------------------------------------------
if isempty(h2d)
    h1 = @(s) 1./(K*abs(PI.k(s)));
    h_minmax = @(s) min(max(h1(s),h_min),h_max);
    [h_gradlimit,H] = GradientLimiting(PI.p(:),h_minmax,g);
    h = h_gradlimit;
    
else % Use pre-computed target mesh size
    h2 = h2d(PI.x(PI.p),PI.y(PI.p));
    h = @(s) interp1(PI.p,h2,s);
end

%--------------------------------------------------------------------------
% Add end points to fixed points
%--------------------------------------------------------------------------
fixedPoints1D = sort(unique(fixedPoints1D(:)));

%--------------------------------------------------------------------------
% Apply force equilibrium
%--------------------------------------------------------------------------
p = ForceEquilibrium1D(0,PI.p(end),h,h_min,fixedPoints1D);

%--------------------------------------------------------------------------
% Store 1D mesh in array
%--------------------------------------------------------------------------
Mesh1D.X = PI.x(p);
Mesh1D.Y = PI.y(p);
Mesh1D.K = PI.k(p);
Mesh1D.p = p;
Mesh1D.h = h;

