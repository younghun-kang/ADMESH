function [Vx,Vy,D,XY1,XY2] = VectorDistanceTransform(XYb,XYq)

Xq = XYq(:,1);
Yq = XYq(:,2);

xyb = cell2mat(XYb(:));
Xb = xyb(:,1);
Yb = xyb(:,2);

IdStart{1} = 1;
IdEnd{1} = length(XYb{1});
for i = 2 : length(XYb)
    IdStart{i} = IdEnd{i-1} + 1;
    IdEnd{i} = IdEnd{i-1} + length(XYb{i});
end

for i = 1 : length(XYb)
    iStart = IdStart{i};
    iEnd = IdEnd{i};
    if isequal(Xb(iStart),Xb(iEnd)) && isequal(Yb(iStart),Yb(iEnd))
        Type{i} = 'Periodic';
    else
        Type{i} = 'nonPeriodic';
    end
end

% Find vector field V1 (from query points to closest boundary node)
id2B = knnsearch([Xb,Yb],[Xq(:),Yq(:)],'k',1);

% Set node id for plus and minus vector
id2Bp = id2B+1;
id2Bm = id2B-1;
    
% Modify node ids of end points based on the line type
for i = 1 : length(Type)
    if strcmpi(Type{i},'nonPeriodic')
        id = id2B == IdStart{i};
        id2Bm(id) = IdStart{i}+1;
        
        id = id2B == IdEnd{i};
        id2Bp(id) = IdEnd{i} - 1;
               
    elseif strcmpi(Type{i},'Periodic')       
        id = id2B == IdStart{i};
        id2Bm(id) = IdEnd{i};
        
        id = id2B == IdEnd{i};
        id2Bp(id) = IdStart{i};
                                      
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

XY1 = [Xq,Yq] + V1;
XY2 = [Xq,Yq] + V1 + V2;

% Restore to original array size
Vx = V(:,1); Vy = V(:,2);
Vx = reshape(Vx,size(Xq,1),size(Xq,2));
Vy = reshape(Vy,size(Xq,1),size(Xq,2));

D = sqrt(Vx.^2 + Vy.^2);


%% ========================================================================
% Debugging block
%==========================================================================
% figure; hold on; axis equal;
% plot(Xb,Yb,'-*k');
% quiver(Xq,Yq,Vx,Vy,0);
% 
% figure; axis equal;
% surf(Xq,Yq,D,'edgecolor','none'); view(2); colorbar; mycmap('redblue');












