addpath(genpath('admesh-lib-part'));
addpath('MeshGeneration1D_v2');

load('test.mat');
temp = struct2table(PTS.Constraints);
Points = temp.xy;
meter2deg = km2deg(1e-3);
h_min          = 10*meter2deg;
h_max          = 100*meter2deg;
delta = h_min/4;

for i = 1 : length(Points)
    x = Points{i}(:,1);
    y = Points{i}(:,2);
    
    [PI(i),CurvLength,p] = ComputeCurveCurvature(x,y);
end

i = 89;
s = linspace(0,1,1e3)';
x = PI(i).x(s);
y = PI(i).y(s);

figure; plot(x,y); axis equal; ax = axis; close; 
% ax = [-93.8503  -93.8497   30.0180   30.0183];
% ax = [-93.8735  -93.8474   30.0006   30.0212];
ax = [-93.8709  -93.8494   30.0012   30.0375];
xl = ax(1); xr = ax(2); yl = ax(3); yr = ax(4);
% xl = min(x); xr = max(x); yl = ax(3); yr = ax(4);

external_boundary = [xl, yl; xr, yl; xr, yr; xl, yr; xl, yl];

PTS = [];
PTS.Poly(1).x = external_boundary(:,1);
PTS.Poly(1).y = external_boundary(:,2);
PTS.Poly(2).x = x;
PTS.Poly(2).y = y;

N = 1e3;
X = linspace(ax(1),ax(2),N);
Y = linspace(ax(3),ax(4),N);
[X,Y] = meshgrid(X,Y);

[D,gradD] = SignedDistanceFunction_v2(PTS,X,Y,delta,h_max);
[Vx,Vy] = Compute8SSED(X,Y,x,y,delta);
[Vx,Vy] = Compute8SSED_v2(x,y,X,Y);

XYb = {external_boundary};
I = [89 90];
for i = 1 : length(PI)
XYb{i+1} = [PI((i)).x(s), PI((i)).y(s)];
end
[Vx,Vy,D] = Compute8SSED_v3(XYb,X,Y,delta);

id = abs(D) < h_min;
figure; hold on; axis equal;
plot(x, y,'k');
quiver(X(id),Y(id),Vx(id),Vy(id),0);

figure; surf(X,Y,D,'EdgeColor','None'); view(2); colorbar; mycmap('bluered');

figure; hold on; axis equal;
for i = 1 : length(PI)
    plot(PI(i).x(s), PI(i).y(s));
end


%%
load('test.mat');
temp = struct2table(PTS.Constraints);
Points = temp.xy;

h_min = Setting1D.h_min;

iMag = 0;
Mag = [];
for i = 1 : length(Points)
    x = Points{i}(:,1);
    y = Points{i}(:,2);
    
    [PI,CurvLength,p] = ComputeCurveCurvature(x,y);
    
    s = linspace(0,1,1e3)';
    x = PI.x(s);
    y = PI.y(s);
    
    nx = -PI.dy(s);
    ny = PI.dx(s);
    nm = sqrt(nx.^2 + ny.^2);
    nx = nx./nm; % normalize
    ny = ny./nm;

    
    
    t = linspace(0,1,1e2)*h_min;
    
    x1 = x + nx*t;
    y1 = y + ny*t;
    
    nx1 = -nx*t;
    ny1 = -ny*t;
        
    % Do for opposite direction
    
    nx = -nx;
    ny = -ny;
    
    x2 = x + nx*t;
    y2 = y + ny*t;
    
    nx2 = -nx*t;
    ny2 = -ny*t;
    
    iMag = iMag + 1;
    Mag.x{iMag,1} = [x1;x2];
    Mag.y{iMag,1} = [y1;y2];
    
    Mag.nx{iMag,1} = [nx1;nx2];
    Mag.ny{iMag,1} = [ny1;ny2];
end

figure; hold on; axis equal;
for i = 1 : length(Mag.x)
    x = Mag.x{i};
    y = Mag.y{i};
    nx = Mag.nx{i};
    ny = Mag.ny{i};
    
    plot(x,y);
    quiver(x,y,nx,ny,0);
    
end

%%
ax = [-93.8503  -93.8497   30.0180   30.0183];
ax = [-93.8735  -93.8474   30.0006   30.0212];
xl = ax(1); xr = ax(2); yl = ax(3); yr = ax(4);

external_boundary = [xl, yl; xr, yl; xr, yr; xl, yr; xl, yl];

PTS = [];
PTS.Poly(1).x = external_boundary(:,1);
PTS.Poly(1).y = external_boundary(:,2);
PTS.Poly(2).x = PI.x(s);
PTS.Poly(2).y = PI.y(s);
delta = h_min;
h_max = Setting1D.h_max;

N = 1e3;
X = linspace(ax(1),ax(2),N);
Y = linspace(ax(3),ax(4),N);
[X,Y] = meshgrid(X,Y);

%**************************************************************************
% IDEA: compute 8ssed at every time step (only for the points near channels)
%**************************************************************************
[D,gradD] = SignedDistanceFunction_v2(PTS,X,Y,delta,h_max);
id = D < 0 & abs(D) < h_min;
[Vx,Vy] = Compute8SSED(X(id),Y(id),PI.x(s),PI.y(s),dx);
figure; hold on; axis equal;
plot(PI.x(s), PI.y(s));
plot(X(id)+Vx,Y(id)+Vy,'.');
% quiver(X(id),Y(id),Vx,Vy,0);

figure; surf(X,Y,D,'edgecolor','none','facecolor','interp'); view(2); axis equal; colorbar;
hold on; myPlotOnTop(PI.x(s),PI.y(s),'k');
mycmap('bluered');

id = D < 0 & abs(D) < h_min;
figure; hold on; axis equal;
plot(PI.x(s), PI.y(s));
quiver(X(id),Y(id),Vx(id),Vy(id),0);

Magx = @(x,y) griddata(X,Y,Vx,x,y);
Magy = @(x,y) griddata(X,Y,Vy,x,y);

N = 1e4;
X = linspace(ax(1),ax(2),N);
Y = linspace(ax(3),ax(4),N);
[X,Y] = meshgrid(X,Y);

[D,gradD] = SignedDistanceFunction_v2(PTS,X,Y,delta,h_max);
id = D > 0 & abs(D) < h_min;
figure; hold on; axis equal;
plot(PI.x(s), PI.y(s),'k');
quiver(X(id),Y(id),Magx(X(id),Y(id)),Magy(X(id),Y(id)),0);

%%
tx1 = PI.dx;
ty1 = PI.dy;
tm = @(s) sqrt(tx1(s).^2 + ty1(s).^2);
nx = @(s) -ty1(s)./tm(s); % normalize
ny = @(s) tx1(s)./tm(s);
nx_h = @(s) nx(s)*h_min;
ny_h = @(s) ny(s)*h_min;

s = linspace(0,1,1e3);
figure; hold on; axis equal;
plot(PI.x(s),PI.y(s));
plot(PI.x(s) + nx_h(s), PI.y(s) + ny_h(s));
% hold on; quiver(PI.x(s) + nx_h(s), PI.y(s) + ny_h(s), -nx_h(s), -ny_h(s),0);

PTS = [];
PTS.Poly(1).x = PI.x(s);
PTS.Poly(1).y = PI.y(s);
PTS.Poly(2).x = PI.x(s) + nx_h(s);
PTS.Poly(2).y = PI.y(s) + ny_h(s);

figure; plot(PI.x(s) + PI.y(s));
hold on; plot(PI.x(s) + PI.y(s) + nx_h(s) + ny_h(s));










