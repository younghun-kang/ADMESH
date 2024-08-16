function MA_new = ConstructMainStream_v3(MA,UIFigure)

msg = 'Constructing main stream...';
progdlg = uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg);

%% ========================================================================
% Construct main stream
%==========================================================================
Size = MA.Size;
BranchNodes = MA.BranchNodes;
BranchNodes1 = MA.BranchNodes;
[temp1,temp2] = cellfun(@(x) ind2sub(MA.Size,x),BranchNodes,'UniformOutput',0);
BranchXY = cellfun(@(x,y) [x(:),y(:)],temp1,temp2,'UniformOutput',0);

BranchNodeEnds = cellfun(@(x) [x(1) x(end)],MA.BranchNodes,'UniformOutput',0);
BranchNodeEnds = cell2mat(BranchNodeEnds(:));

% for i = 1 : length(BranchConnectivity)
%     ID = BranchNodes{i};
%     [I,J] = ind2sub(Size,ID);
%     BranchLength(i) = sum(sqrt( (I(1:end-1) - I(2:end)).^2 + (J(1:end-1) - J(2:end)).^2));
% end


NodeEndList = unique(BranchNodeEnds(:));
MergeID = zeros(length(NodeEndList),2);
for k = 1 : length(NodeEndList)
%     if k < 1579
%         continue;
%     end
    I1 = find(BranchNodeEnds(:,1) == NodeEndList(k));
    I2 = find(BranchNodeEnds(:,2) == NodeEndList(k));
    I = [I1; I2];
    if length(I) == 1
        continue;
    elseif length(I) == 2
%         MergeID(k,:) = I;
%         continue;
    end
    if ~isempty(I1) || ~isempty(I2)
        x = cell(length(I1)+length(I2),1);
        y = cell(length(I1)+length(I2),1);
        for i = 1 : length(I1)
            x{i} = BranchXY{I1(i)}(:,1);
            y{i} = BranchXY{I1(i)}(:,2);
        end
        for j = 1 : length(I2)
            x{length(I1) + j} = flipud(BranchXY{I2(j)}(:,1));
            y{length(I1) + j} = flipud(BranchXY{I2(j)}(:,2));
        end
    end

    x1 = x;
    y1 = y;
%     x1 = cellfun(@(x) x(1:min(10,length(x))),x1,'UniformOutput',false);
%     y1 = cellfun(@(x) x(1:min(10,length(x))),y1,'UniformOutput',false);
    
    SmoothingRMSE = km2deg(.5e-3);
    K = inf(length(x1));
%     figure; hold on;
    for i = 1 : length(x1)
        for j = i+1 : length(x1)
            x2 = [flip(x1{j}); x1{i}(2:end)];
            y2 = [flip(y1{j}); y1{i}(2:end)];
           
            FixedPoints = [x2([1,end]),y2([1,end]); x1{i}(1), y1{i}(1)];
            PI = ComputePathCurvature([x2,y2],FixedPoints,SmoothingRMSE);
            id = length(x1{j});
            K(i,j) = abs(PI.k(PI.p(id)));
%             plot(PI.x1(PI.p),PI.y1(PI.p));
%             myPlot3(PI.sx(PI.p),PI.sy(PI.p),PI.k(PI.p));
        end
    end
    [i,j] = find(K == min(K(:)),1);
    MergeID(k,:) = I([i,j]);
    
%     SmoothingRMSE = km2deg(.5e-3);
%     dx = zeros(length(x),1);
%     dy = zeros(length(x),1);
% %     figure; hold on;
%     for i = 1 : length(x1)
%         FixedPoints = [x1{i}([1,end]),y1{i}([1,end])];
%         PI = ComputePathCurvature([x1{i},y1{i}],FixedPoints,SmoothingRMSE);
%         if isempty(PI)
%             dx(i) = x1{i}(2) - x1{i}(1);
%             dy(i) = y1{i}(2) - y1{i}(1);
%         else
%             dx(i) = PI.sx(PI.p(2)) - PI.sx(PI.p(1));
%             dy(i) = PI.sy(PI.p(2)) - PI.sy(PI.p(1));
%         end
% %         plot(x1{i},y1{i});
% %         plot(PI.sx(PI.p),PI.sy(PI.p));
%     end
%     s = sqrt(dx.^2 + dy.^2);
%     dx = dx./s;
%     dy = dy./s;
%     
%     DotProd = dx*dx' + dy*dy';
%     [i,j] = find(DotProd == min(DotProd(:)),1);
%     if min(DotProd(:)) < -0.5
%         MergeID(k,:) = I([i,j]);
%     end
    
% %     SmoothingParam = 1;
%     SmoothingRMSE = km2deg(1e-4);
%     dx = zeros(length(x),1);
%     dy = zeros(length(x),1);
% %     figure; axis equal;
%     for i = 1 : length(x)
%         x1 = x{i};
%         y1 = y{i};
% %         x1 = x1(1:min(10,length(x1)));
% %         y1 = y1(1:min(10,length(y1)));
%         t = 1 : length(x1);
% %         x1 = csaps(t,x1,SmoothingParam,t);
% %         y1 = csaps(t,y1,SmoothingParam,t);
%         x1 = csapsRMSE(t,x1,SmoothingRMSE,t);
%         y1 = csapsRMSE(t,y1,SmoothingRMSE,t);
%         dx(i) = x1(1) - x1(2);
%         dy(i) = y1(1) - y1(2);
% %         hold on; plot(x1,y1);
% %         quiver(x1(2),y1(2),dx(i),dy(i));
%     end
% %     axis([-94.1063  -94.1063   30.1394   30.1394]);
%     
%     s = sqrt(dx.^2 + dy.^2);
%     dx = dx./s;
%     dy = dy./s;
%     
%     temp = dx*dx' + dy*dy';
%     [i,j] = find(temp == min(temp(:)),1);


    
    progdlg.Value = k/length(NodeEndList);
end
close(progdlg);

MergeID1 = MergeID;
BranchNodes1 = BranchNodes;
for k = 1 : size(MergeID1,1)
    I1 = MergeID1(k,1);
    I2 = MergeID1(k,2);
    if I1 == I2
        continue;
    end
    
    nodes1 = BranchNodes1{I1};
    nodes2 = BranchNodes1{I2};
    nodes1 = nodes1(:);
    nodes2 = nodes2(:);
    
    if nodes1(end) == nodes2(1)
        nodes = [nodes1; nodes2(2:end)];
    elseif nodes1(end) == nodes2(end)
        nodes = [nodes1; flipud(nodes2(1:end-1))];
    elseif nodes1(1) == nodes2(1)
        nodes = [flipud(nodes1); nodes2(2:end)];
    elseif nodes1(1) == nodes2(end)
        nodes = [flipud(nodes1); flipud(nodes2(1:end-1))];
    else
        % This case happens when a loop
        continue;
    end
    
    BranchNodes1{I1} = nodes;
    BranchNodes1{I2} = [];
    
    J1 = MergeID1 == I2;
    MergeID1(J1) = I1;
end

I = cellfun(@(x) isempty(x),BranchNodes1);
BranchNodes1(I) = [];
BranchNodes1 = BranchNodes1(:);


MA_new.Size = MA.Size;
MA_new.BranchNodes = BranchNodes1;
% MA_new.BranchLength = MainStreamLength;

return;

figure; hold on; axis equal;
for i = 1 : length(BranchXY1)
%     plot(BranchXY{i}(:,1),BranchXY{i}(:,2));
    plot(BranchXY1{i}(:,1),BranchXY1{i}(:,2));
end


MainStreamBranch = {};
MainStreamLength = [];
wbar = waitbar(0);
for i = 1 : length(BranchConnectivity)
    
    % Case 1: Isolated branch
    if isempty(BranchConnectivity{i})
        MainStreamBranch{end+1} = i;
        MainStreamLength(end+1) = BranchLength(i);
        continue;
    end
    
    x1 = BranchXY{i}(:,1);
    y1 = BranchXY{i}(:,2);
    t1 = 1 : length(x1);
    
    x1 = csaps(t1,x1,1e-2,t1);
    y1 = csaps(t1,y1,1e-2,t1);
    
    for j = 1 : length(BranchConnectivity{i})
        k = BranchConnectivity{i}(j);
        x2 = BranchXY{k}(:,1);
        y2 = BranchXY{k}(:,2);
        t2 = 1 : length(x2);
    
        x2 = csaps(t2,x2,1e-2,t2);
        y2 = csaps(t2,y2,1e-2,t2);
        
    end
    
    % Case 2: Not a level-1 branch
    if ~all(BranchConnectivity{i} > 0) && ~all(BranchConnectivity{i} < 0)
        continue;
    end
    
    % Case 3: Level 1 branch
    BranchChunk = i;
    while 1
        ChunkAdd = [];
        for j = 1 : length(BranchChunk)
            k = BranchChunk(j);
            ChunkAdd = [ChunkAdd; setdiff(abs(BranchConnectivity{k}),BranchChunk)];
        end
        if isempty(ChunkAdd)
            break;
        end
        BranchChunk = [BranchChunk; ChunkAdd];
    end
    
    for j = 1 : length(BranchChunk)
        k = BranchChunk(j);
        if all(BranchConnectivity{k} > 0) || all(BranchConnectivity{k} < 0)
            BranchChunk(j) = -BranchChunk(j);
        end
    end
    BranchChunk = sort(BranchChunk);
    
    
    while 1
        Path2 = [];
        for j = 1 : nnz(BranchChunk < 0)
            k = -BranchChunk(j);
            Path = {k};
            while 1
                PathNew = [];
                for kk = 1 : length(Path)
                    PathAdd = setdiff(abs(BranchConnectivity{Path{kk}(end)}),unique(vertcat(Path{:})));
                    PathAdd = setdiff(PathAdd,vertcat(MainStreamBranch{:}));
                    %             PathAdd = setdiff(PathAdd,unique([Path{:}]));
                    if isempty(PathAdd)
                        PathNew{end+1} = Path{kk}(:);
                    else
                        for kkk = 1 : length(PathAdd)
                            PathNew{end+1} = [Path{kk}(:); PathAdd(kkk)];
                        end
                    end
                end
                if isequal(size(Path),size(PathNew)) && all(cellfun(@(x,y) all(ismember(x,y)),Path,PathNew))
                    break;
                end
                Path = PathNew;
            end
            Path2 = [Path2, Path];
        end
        
        for j = 1 : length(Path2)
            Path2{j} = setdiff(Path2{j},vertcat(MainStreamBranch{:}),'stable');
        end
        
        PathLength = zeros(length(Path2),1);
        for j = 1 : length(Path2)
            PathLength(j) = sum(BranchLength(Path2{j}));
        end
        [maxlength,imax] = max(PathLength);
        
        if isempty(Path2{imax})
            break;
        end
        MainStreamBranch{end+1} = Path2{imax};
        MainStreamLength(end+1) = maxlength;
    end
    waitbar(i/length(BranchConnectivity),wbar);
end
delete(wbar);

%--------------------------------------------------------------------------
% Handle not included branches... especially short bridges
%--------------------------------------------------------------------------
id = find(~ismember(1:length(BranchConnectivity),vertcat(MainStreamBranch{:})));
temp_Length = BranchLength(id);
id = mat2cell(id(:),ones(length(id),1),1);

MainStreamBranch = [MainStreamBranch(:); id(:)];
MainStreamLength = [MainStreamLength(:); temp_Length(:)];


BranchNodes_new = {};
MainStreamNodes = cell(length(MainStreamBranch),1);

wbar = waitbar(0);
for i = 1 : length(MainStreamBranch)
    if length(MainStreamBranch{i}) == 1
        k = MainStreamBranch{i};
        MainStreamNodes{i} = BranchNodes{k};
        JointConnectivity_new(i,:) = [BranchNodes{k}(1),BranchNodes{k}(end)];
        continue;
    end
    
    k = MainStreamBranch{i}(1);
    k1 = MainStreamBranch{i}(2);
    if any(BranchConnectivity{k} == k1)
        iMainStreamNodes = (BranchNodes{k});
    elseif any(BranchConnectivity{k} == -k1)
        iMainStreamNodes = flip(BranchNodes{k});
    end
        
    for j = 2 : length(MainStreamBranch{i})
        k1 = MainStreamBranch{i}(j-1);
        k = MainStreamBranch{i}(j);
        if any(BranchConnectivity{k} == -k1)
            iMainStreamNodes = [iMainStreamNodes; (BranchNodes{k})];
        elseif any(BranchConnectivity{k} == k1)
            iMainStreamNodes = [iMainStreamNodes; flip(BranchNodes{k})];
        end
    end
    
    iMainStreamNodes = unique(iMainStreamNodes,'rows','stable');
    MainStreamNodes{i} = iMainStreamNodes;
    JointConnectivity_new(i,:) = [iMainStreamNodes(1),iMainStreamNodes(end)];
    waitbar(i/length(MainStreamBranch),wbar);
end
delete(wbar);

MA_new.Size = MA.Size;
MA_new.nBranch = length(MainStreamNodes);
MA_new.BranchNodes = MainStreamNodes;
MA_new.JointConnectivity = JointConnectivity_new;
MA_new.BranchLength = MainStreamLength;
% MA_new.XY = XY_new;

return;


JointConnectivity1 = MA.JointConnectivity;
XY = MA.XY;
id_channels = [];
ICONN = [];

BranchNodes_new = BranchNodes;
JointConnectivity_new = [];
XY_new = XY;
max_length = zeros(length(BranchNodes_new),1);
i = 0;
while 1
    if i >= size(JointConnectivity,1)
        break;
    end
    i = i + 1;
    
    % skip if both ends are not a tip
    if nnz(JointConnectivity == JointConnectivity(i,1)) ~= 1 && nnz(JointConnectivity == JointConnectivity(i,2)) ~= 1 
        continue;
    end
    
    num_channels = 0;
    id_channels = [];
    
    % Make the array "JointConnectivity" stores id in left-to-right order
    if nnz(JointConnectivity == JointConnectivity(i,2)) == 1
        JointConnectivity(i,:) = JointConnectivity(i,[2 1]);
%         XY{i} = flip(XY{i});
    end
    
    if nnz(JointConnectivity == JointConnectivity(i,2)) == 1 % skip if tip-to-tip
        continue;
    end
    
    num_channels = num_channels + 1;
    id_channels{num_channels} = JointConnectivity(i,:);
    
    JointConnectivity(i,:) = nan;
%     BranchNodes(i) = nan;
    k = 0;
    while 1
        k = k + 1;
        if k > length(id_channels)
            break;
        end
        ic = id_channels{k};
        
        j1 = find(JointConnectivity(:,1) == ic(end));
        j2 = find(JointConnectivity(:,2) == ic(end));
        JointConnectivity(j2,:) = JointConnectivity(j2,[2 1]);
%         for jj = 1 : length(j2)
%             XY{j2(jj)} = flip(XY{j2(jj)});
%         end
        j = vertcat(j1,j2);
        
        for jj = 1 : length(j)
            jjj = j(jj);
            if JointConnectivity(jjj,2) > 0
                id_channels{end+1} = [ic,JointConnectivity(jjj,2)];
            else
                id_channels{end+1} = [ic,JointConnectivity(jjj,2),nan];
            end
        end
        JointConnectivity(j,:) = nan;
%         BranchNodes(j) = nan;
    end
    
    % Remove channel connections which are entirely included in others
    id_channels = flip(id_channels);
    iremove = [];
    for j = 1 : length(id_channels)
        for k = j+1 : length(id_channels)
            if all(ismember(id_channels{k},id_channels{j}))
                id_channels{k} = [];
                iremove = [iremove,k];
            end
        end
    end
    id_channels(iremove) = [];

    % Compute length of channels
    LENGTH = zeros(length(id_channels),1);
    ICONN = cell(length(id_channels),1);
    for j = 1 : length(id_channels)
        id = id_channels{j};
        id(isnan(id)) = [];
        ICONN{j} = [];
        for k = 1 : (length(id)-1)
            kk = find(ismember(JointConnectivity1,id([k k+1]),'rows'));
            if isempty(kk)
                kk = -find(ismember(JointConnectivity1,id([k+1 k]),'rows'));
            end
            ICONN{j} = [ICONN{j}, kk];
            id1 = BranchNodes1{abs(kk)};
            [I,J] = ind2sub(Size,id1);
            dist = sum(sqrt( (I(1:end-1) - I(2:end)).^2 + (J(1:end-1) - J(2:end)).^2));
            LENGTH(j) = LENGTH(j) + dist;
        end
    end
    
    [~,id_main] = max(LENGTH);
    XY_add = [];
    id_new = [];
    ICONN_main = ICONN{id_main};
    for j = 1 : length(ICONN_main)
        if ICONN_main(j) > 0
            id_add = BranchNodes1{ICONN_main(j)};
            XY_add = vertcat(XY_add,XY{ICONN_main(j)});
        else
            id_add = flip(BranchNodes1{-ICONN_main(j)});
            XY_add = vertcat(XY_add,flip(XY{-ICONN_main(j)}));
        end
        id_new = vertcat(id_new,id_add);
        BranchNodes_new{abs(ICONN_main(j))} = [];
        
        
        
    end
    if ~isempty(id_new)
        BranchNodes_new{end+1} = id_new;
        XY_new{end+1} = XY_add;
        max_length(end+1) = max(LENGTH);
    end
    
%     BranchNodes_new = vertcat(BranchNodes{ICONN{id_main}});
    
%     BranchNodes{ICONN{id_main}(1)} = BranchNodes_new;
%     JointConnectivity1(ICONN{id_main}(1),:) = [BranchNodes_new(1), BranchNodes_new(end)];
    
%     JointConnectivity1(ICONN{id_main}(2:end),:) = 0;
    
    0;
    
    
end

id_remove = [];
for i = 1 : length(BranchNodes_new)
    BranchNodes_new{i} = unique(BranchNodes_new{i},'stable');
    XY_new{i} = unique(XY_new{i},'rows','stable');
    if isempty(BranchNodes_new{i})
        id_remove = [id_remove,i];
    end
end
BranchNodes_new(id_remove) = [];
XY_new(id_remove) = [];

for i = 1 : length(BranchNodes_new)
    JointConnectivity_new(i,:) = [BranchNodes_new{i}(1),BranchNodes_new{i}(end)];
end

% XY_new = cell(length(BranchNodes_new),1);
% for i = 1 : length(BranchNodes_new)
%     XY_new{i} = [X(BranchNodes_new{i}),Y(BranchNodes_new{i})];
% end

MA_new.Size = MA.Size;
MA_new.nBranch = length(BranchNodes_new);
MA_new.BranchNodes = BranchNodes_new;
MA_new.JointConnectivity = JointConnectivity_new;
MA_new.XY = XY_new;


