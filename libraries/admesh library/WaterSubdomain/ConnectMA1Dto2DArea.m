function MA_connected = ConnectMA1Dto2DArea(MA,L2D_new_boundary,AreaInfo,X,Y,dx)

%% ========================================================================
% Make MA1D is connected to 2D area
%==========================================================================
XY_L2D_boundary = vertcat(L2D_new_boundary{:});
L1D_Filled = AreaInfo.L1D;

nBranch = MA.nBranch;
BranchNodes = MA.BranchNodes;
JointConnectivity = MA.JointConnectivity;
XY = cell(length(BranchNodes),1);
Length = zeros(length(BranchNodes),1);

if isempty(L2D_new_boundary)
    for i = 1 : nBranch
        id = BranchNodes{i};
        XY{i} = [X(id),Y(id)];
    end
    MA_connected = MA;
    MA_connected.BranchXY = XY;
    return;
end

for i = 1 : nBranch
    id = BranchNodes{i};
    k = L1D_Filled(id) > 0;
    id = id(k);
    if isempty(id)
        continue;
    end
    x = X(id);
    y = Y(id);
    
    dist1 = sqrt((XY_L2D_boundary(:,1) - x(1)).^2 + (XY_L2D_boundary(:,2) - y(1)).^2);
    dist2 = sqrt((XY_L2D_boundary(:,1) - x(end)).^2 + (XY_L2D_boundary(:,2) - y(end)).^2);
    
    if any(dist1 < (2)*dx)
        [~,iend] = min(dist1);
%         MA_connected.ID{i} = vertcat(iend,MA_connected.ID{i});
        x = vertcat(XY_L2D_boundary(iend,1),x);
        y = vertcat(XY_L2D_boundary(iend,2),y);
        JointConnectivity(i,1) = -JointConnectivity(i,1);
    end
    
    if any(dist2 < (2)*dx)
        [~,iend] = min(dist2);
%         MA_connected.ID{i} = vertcat(MA_connected.ID{i},iend);
         x = vertcat(x,XY_L2D_boundary(iend,1));
         y = vertcat(y,XY_L2D_boundary(iend,2));
         JointConnectivity(i,2) = -JointConnectivity(i,2);
    end
    XY{i} = [x,y];
    Length(i) = sum( sqrt((x(1:end-1) - x(2:end)).^2 + (y(1:end-1) - y(2:end)).^2));
end
[~,id] = unique(JointConnectivity,'rows');

MA_connected.nBranch = length(id);
MA_connected.BranchNodes = BranchNodes(id);
MA_connected.JointConnectivity = JointConnectivity(id,:);
MA_connected.BranchXY = XY(id);
MA_connected.BranchLength = Length(id);

