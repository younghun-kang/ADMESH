function MA_smoothed = SmoothingMA(MA,SmoothingParam)

% MA = MA_local;
BranchNodes = MA.BranchNodes;
JointConnectivity = MA.JointConnectivity;

BranchXY = MA.BranchXY;
BranchXY_smooth = cell(MA.nBranch,1);

for i = 1 : MA.nBranch
    xy = BranchXY{i};
    
    id = false(size(xy,1),1);
%     for j = 1 : size(xy,1)
%         id(j) = nnz(ismember(vertcat(BranchXY{:}),xy(j,:),'rows')) > 1;
%     end
    
    [~,temp] = ismember(vertcat(BranchXY{:}),xy,'rows');
    for j = 1 : size(xy,1)
        id(j) = nnz(temp == j) > 1;
    end
    
    if JointConnectivity(i,1) < 0
        id(1) = true;
    end
    if JointConnectivity(i,2) < 0
        id(end) = true;
    end
    FixedPointsXY = xy(id,:);
%     id = ismember(xy,vertcat(XY{:}),'rows');
    
%     id = [1, size(xy,1)];
%     k = [Conn(i,1) < 0 || nnz(Conn == Conn(i,1)) > 1;
%         Conn(i,2) < 0 || nnz(Conn == Conn(2,1)) > 1];

%     FixedPointsXY = xy(id(k),:);
%     FixedPointsXY = [];
    
    id_fixedPoints = [];
    if ~isempty(FixedPointsXY)
        id_fixedPoints = find(ismember(xy,FixedPointsXY,'rows'));
    end
    id_fixedPoints = [1; id_fixedPoints(:); size(xy,1)];
    
    x = xy(:,1); y = xy(:,2);
    
    t = 1 : length(x);
    w = ones(length(x),1);
    w(id_fixedPoints) = 1e8;
%     SmoothingParam = min(0.01,sqrt(1/length(t)));

    Smooth_x = csaps(t,x,SmoothingParam,t,w);
    Smooth_x(id_fixedPoints) = x(id_fixedPoints);
    
    Smooth_y = csaps(t,y,SmoothingParam,t,w);
    Smooth_y(id_fixedPoints) = y(id_fixedPoints);

    BranchXY_smooth{i} = [Smooth_x(:),Smooth_y(:)];
    fprintf('Smoothing MA (%.2f)\n',i/MA.nBranch);
end

MA_smoothed = MA;
MA_smoothed.BranchXY = BranchXY_smooth;



