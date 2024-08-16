function [Vx,Vy,D] = Compute8SSED_v2(Xb,Yb,Xq,Yq)

if isequal(Xb(1),Xb(end)) && isequal(Yb(1),Yb(end))
    Type = 'Periodic';
else
    Type = 'nonPeriodic';
end

Debugging = 0;

% Find vector field V1 (from query points to closest boundary node)
id2B = knnsearch([Xb,Yb],[Xq(:),Yq(:)],'k',1);
V1 = [Xb(id2B) - Xq(:), Yb(id2B) - Yq(:)];

% Set node id for plus and minus vector
id2Bp = id2B+1;
id2Bm = id2B-1;
    
% Modify node ids of end points based on the line type
if strcmpi(Type,'nonPeriodic')
    id = id2B == 1;
    id2Bm(id) = 2;
    
    id = id2B == length(Xb);
    id2Bp(id) = length(Xb) - 1;
    
    % Store end points
%     V(id,:) = V1(id,:);
%     V1(id,:) = []; id2B(id) = [];
    
elseif strcmpi(Type,'Periodic')
    id2Bp = id2B+1;
    id2Bm = id2B-1;
    
    id = id2B == 1;
    id2Bm(id) = length(Xb);
    
    id = id2B == length(Xb);
    id2Bp(id) = 1;
end

Vp = [Xb(id2Bp) - Xb(id2B), Yb(id2Bp) - Yb(id2B)];
Vm = [Xb(id2Bm) - Xb(id2B), Yb(id2Bm) - Yb(id2B)];

% Compute dot products
DotProd = @(v1,v2) sum(v1.*v2,2);

V1Vp = DotProd(V1,Vp);
V1Vm = DotProd(V1,Vm);
VpVp = DotProd(Vp,Vp);
VmVm = DotProd(Vm,Vm);

% Check concave part. Take edge where has lower angle.
id = V1Vp < 0 & V1Vm < 0 & V1Vp./sqrt(VpVp) > V1Vm./sqrt(VmVm);
V2 = Vp;
V2(id,:) = Vm(id,:);

% Check angle between Vp and Vm
id = V1Vp > 0 & V1Vm < 0;
V2(id,:) = Vm(id,:);

% Compute t parameter and V
V1V2 = DotProd(V1,V2);
V2V2 = DotProd(V2,V2);

t = - V1V2./(V2V2);
t(t < 0) = 0; % negative t occurs for convex part

V = V1 + repmat(t,1,2).*V2;

% Restore to original array size
Vx = V(:,1); Vy = V(:,2);
Vx = reshape(Vx,size(Xq,1),size(Xq,2));
Vy = reshape(Vy,size(Xq,1),size(Xq,2));

D = sqrt(Vx.^2 + Vy.^2);
if Debugging == 1
    figure; hold on; axis equal;
    plot(Xb,Yb,'-*k');
    quiver(Xq,Yq,Vx,Vy,0);
end












