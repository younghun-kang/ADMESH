function [Vx,Vy,D] = Compute8SSED_v4(ipgon,Xq,Yq,dx)

I = find(isnan(ipgon.Vertices(:,1)));
I = [0; I(:); size(ipgon.Vertices,1)+1];
XYb = [];
for i = 1 : length(I)-1
    i1 = I(i) + 1;
    i2 = I(i+1) - 1;
    XYb{i} = ipgon.Vertices([i1:i2, i1],:);
end

for i = 1 : length(XYb)
    PTS.Poly(i).x = XYb{i}(:,1);
    PTS.Poly(i).y = XYb{i}(:,2);
end
% IN = PointsInDomain2(Xq,Yq,PTS);

for i = 1 : length(XYb)
    x = XYb{i}(:,1);
    y = XYb{i}(:,2);
    t = 1 : length(x);
    
    dist = sqrt((x(1:end-1)-x(2:end)).^2 + (y(1:end-1)-y(2:end)).^2);
    if isempty(dx)
        dx = mean(mean(diff(Xq,1,2)));
    end
    
%     if all(dist < dx); continue; end
    
    t2 = 1;
    for j = 1 : length(x)-1
        n = max(2,round(dist(j)/dx));
        t2 = [t2, linspace(t2(end),t2(end)+1,n)];
    end
    t2 = unique(t2);
    t2 = t2(:);
    
    XYb{i} = [interp1(t,x,t2), interp1(t,y,t2)];

end

xyb = cell2mat(XYb(:));
Xb = xyb(:,1);
Yb = xyb(:,2);

Type = {1, length(XYb{1})};
for i = 2 : length(XYb)
    Type{i,1} = Type{i-1,2} + 1;
    Type{i,2} = Type{i-1,2} + length(XYb{i});
end

for i = 1 : length(XYb)
    iStart = Type{i,1};
    iEnd = Type{i,2};
    if isequal(Xb(iStart),Xb(iEnd)) && isequal(Yb(iStart),Yb(iEnd))
        Type{i,3} = 'Periodic';
    else
        Type{i,3} = 'nonPeriodic';
    end
end

% Find vector field V1 (from query points to closest boundary node)
id2B = knnsearch([Xb,Yb],[Xq(:),Yq(:)],'k',1);

% Set node id for plus and minus vector
id2Bp = id2B+1;
id2Bm = id2B-1;
    
% Modify node ids of end points based on the line type
for i = 1 : size(Type,1)
    if strcmpi(Type{i,3},'nonPeriodic')
        id = id2B == Type{i,1};
        id2Bm(id) = Type{i,1}+1;
        
        id = id2B == Type{i,2};
        id2Bp(id) = Type{i,2} - 1;
               
    elseif strcmpi(Type{i,3},'Periodic')       
        id = id2B == Type{i,1};
        id2Bm(id) = Type{i,2};
        
        id = id2B == Type{i,2};
        id2Bp(id) = Type{i,1};
    end
end

V1 = [Xb(id2B) - Xq(:), Yb(id2B) - Yq(:)];
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
t(isnan(t)) = 0; % NaN occurs when duplicated points exist

V = V1 + repmat(t,1,2).*V2;

% Restore to original array size
Vx = V(:,1); Vy = V(:,2);
Vx = reshape(Vx,size(Xq,1),size(Xq,2));
Vy = reshape(Vy,size(Xq,1),size(Xq,2));

D = sqrt(Vx.^2 + Vy.^2);
% D(IN) = -D(IN);


%% ========================================================================
% Debugging block
%==========================================================================
% figure; hold on; axis equal;
% plot(Xb,Yb,'-*k');
% quiver(Xq,Yq,Vx,Vy,0);
% 
% figure; axis equal;
% surf(Xq,Yq,D,'edgecolor','none'); view(2); colorbar; mycmap('redblue');












