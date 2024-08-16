function MA_pruned = SkeletonPruning_Choi_level1branch(Xq,Yq,Vx,Vy,D,MA,rho_pruning,dtheta)
Debugging = 0;

BranchNodes   = MA.BranchNodes;
JointConnectivity = MA.JointConnectivity;
if isfield(MA,'PruningFlag')
    PruningFlag = MA.PruningFlag;
end
% XY   = MA.XY;

iBranch_pruned = [];
ID_pruned = [];

N = size(Xq,1);
M = size(Xq,2);
dx = mean(mean(diff(Xq,1,2)));
dy = mean(mean(diff(Yq,1,1)));

while 1
    nBranch_pruned = 0;
    for i = 1 : length(BranchNodes)
        num_nghb(1) = nnz(JointConnectivity == JointConnectivity(i,1));
        num_nghb(2) = nnz(JointConnectivity == JointConnectivity(i,2));
        id_branch_tip = find(num_nghb == 1);
        
        if isfield(MA,'PruningFlag') && PruningFlag(i) == 0
            id_branch_tip = [];
        end
        
        if isempty(id_branch_tip) % Keep all points if it's not a level-1 branch
            id_new = BranchNodes{i};
%             XY_new = XY{i};
        elseif length(id_branch_tip) > 1 % Remove disconnected branches
            id_new = [];
        elseif length(id_branch_tip) == 1
            id_new = BranchNodes{i};
%             XY_new = XY{i};
            
            k = id_new;
            [I,J] = ind2sub(size(Xq),k);
            k(I == 1 | I == size(Xq,1) | J == 1 | J == size(Xq,2)) = [];
            
            k_nghb = k + [-N-1, -N, -N+1, -1, 1, N-1, N, N+1];
            
            Dsq = (Vx(k_nghb) - Vx(k)).^2 + (Vy(k_nghb) - Vy(k)).^2;
            Qi_norm = Vx(k_nghb).^2 + Vy(k_nghb).^2;
            Q_norm = Vx(k).^2 + Vy(k).^2;
            
            ConnectivityTest = (Qi_norm - Q_norm) ./ max(Xq(k) - Xq(k_nghb),Yq(k)-Yq(k_nghb));
            
            COSINE = (Vx(k_nghb).*Vx(k) + Vy(k_nghb).*Vy(k))...
                ./(sqrt(Vx(k_nghb).^2 + Vy(k_nghb).^2) .* sqrt(Vx(k).^2 + Vy(k).^2));
            COSINE(isnan(COSINE)) = 1;
            theta = acos(COSINE);
            
            if length(k) == 1
                max_Dsq = max(Dsq);
                max_theta = max(theta);
            else
                max_Dsq = max(Dsq,[],2);
                max_theta = max(theta,[],2);
            end
            
            %----------------------------------------------------------------------
            % Prune if skeleton points are near corner
            %----------------------------------------------------------------------
            id_corner = sqrt(max_Dsq) < rho_pruning & max_theta <= dtheta;
            id_remove = id_corner;
            %----------------------------------------------------------------------
            % Prune if skeleton points are too close to boundary
            %----------------------------------------------------------------------
            id_boundary = abs(D(k)) < sqrt(dx^2 + dy^2) | max_theta < 1e-8;
            id_remove = id_remove | id_boundary;
            
            if id_branch_tip == 1
                j = find(id_remove == false,1,'first');
                id_remove(j:end) = false;
                
            else
                j = find(id_remove == false,1,'last');
                id_remove(1:j) = false;
                
            end
            id_new(id_remove) = [];
%             XY_new(id_remove,:) = [];
        end
        
        if ~isempty(id_new) && length(id_new) > 1
            nBranch_pruned = nBranch_pruned + 1;
            ID_pruned{nBranch_pruned} = id_new;
%             XY_pruned{nBranch_pruned} = XY_new;
        else
            iBranch_pruned = vertcat(iBranch_pruned,i);
        end
    end
    if isequal(BranchNodes,ID_pruned) && isempty(iBranch_pruned)
        break;
    end
    BranchNodes = ID_pruned;
    ID_pruned = [];
%     XY = XY_pruned;
%     XY_pruned = [];
    JointConnectivity(iBranch_pruned,:) = [];
    if isfield(MA,'PruningFlag') 
        PruningFlag(iBranch_pruned) = [];
    end
    
    iBranch_pruned = [];
end

MA_pruned.nBranch = nBranch_pruned;
MA_pruned.Size = MA.Size;
MA_pruned.BranchNodes = BranchNodes;
MA_pruned.JointConnectivity = JointConnectivity;
% MA_pruned.XY = XY;

if Debugging == 1
    figure; hold on;
    set(gcf,'position',[1922 258 766 738]);
    for i = 1 : length(MA.ID)
        plot(Xq(MA.ID{i}),Yq(MA.ID{i}),'linewidth',1.5);
    end
    axis equal;
    
    figure; hold on;
    set(gcf,'position',[2690 258 766 738]);
    for i = 1 : length(ID_pruned)
        plot(Xq(ID_pruned{i}),Yq(ID_pruned{i}),'linewidth',1.5);
    end
    axis equal;
    
    figure; hold on;
    for i = 1 : length(ID_pruned)
        plot(XY{i}(:,1),XY{i}(:,2),'linewidth',1.5);
    end
    axis equal;
end















