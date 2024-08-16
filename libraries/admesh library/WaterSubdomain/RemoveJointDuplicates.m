function MA_new = RemoveJointDuplicates(MA)

Size = MA.Size;

BranchNodes = MA.BranchNodes;
BranchLength = cellfun(@(x) length(x),BranchNodes);
BranchList = find(BranchLength(:) == 2);

% Find branches with length = 2
BranchNodeEnds = cellfun(@(x) [x(1) x(end)],BranchNodes,'UniformOutput',0);
BranchNodeEnds = cell2mat(BranchNodeEnds(:));

BNE = BranchNodeEnds(BranchList,:);

% Loop over branches
for i = 1 : length(BranchList)
    % Repeat for each end point
    for j = 1 : 2
        iBranch = i;
        iJ = BNE(iBranch,j);
        
        PathBI = {[i,j]}; % Path in terms of branch indexes
        PathBNE = {iJ}; % Path in terms of branch node ends
        
        if iJ == 187691581
            0;
        end
        %------------------------------------------------------------------
        % Find all possible paths of joints
        %------------------------------------------------------------------
        while 1
            nPath = length(PathBI);
            PathBI_old = PathBI;
            
            for k = 1 : nPath
                % Retrieve the last joint id
                iBranch = PathBI{1}(end,1);
                jBranch = PathBI{1}(end,2);
                
                % Find branches sharing the new end point
                iJ = BNE(iBranch,jBranch);
                
                [II,JJ] = find(BNE == iJ);
                JJ(II == iBranch) = [];
                II(II == iBranch) = [];
                
                % Switch j index to find end point on the other side
                j1 = JJ == 1;
                j2 = JJ == 2;
                JJ(j1) = 2;
                JJ(j2) = 1;
                
                K = sub2ind(size(BNE),II,JJ);
                iJ = BNE(K);
                
                % Add new joint in the path
                if isempty(II) || all(ismember(II,PathBI{1}(2:end,1)))
                    PathBI{end+1} = PathBI{1};
                    PathBNE{end+1} = PathBNE{1};
                else
                    for kk = 1 : length(II)
                        if ismember(II(kk),PathBI{1}(2:end,1))
                            continue;
                        end
                        PathBI{end+1} = [PathBI{1}; II(kk), JJ(kk)];
                        PathBNE{end+1} = [PathBNE{1}; iJ(kk)];
                    end
                end
                PathBI(1) = [];
                PathBNE(1) = [];
                
                PathBI = PathBI(:);
                PathBNE = PathBNE(:);
            end
            if isequal(PathBI_old,PathBI)
                break;
            end
        end
        
%         %------------------------------------------------------------------
%         % If there is only one possible path, it means 
%         %------------------------------------------------------------------
%         if length(PathBI) == 1
%             continue;
%         end
        
        %------------------------------------------------------------------
        % Find paths with a loop
        %------------------------------------------------------------------
        I = cellfun(@(x) x(1) == x(end),PathBNE);
        PathBNE = PathBNE(I);
        PathBI = PathBI(I);
        
        if isempty(PathBI) % If no loop
            continue;
        end
        
        %------------------------------------------------------------------
        % Remove pathes entirely included in another one
        %------------------------------------------------------------------
        I = cellfun(@(x) length(x),PathBNE);
        [~,I] = sort(I);
        PathBNE = PathBNE(I);
        PathBI = PathBI(I);
        
        for k = 1 : length(PathBNE)
            for kk = k+1:length(PathBNE)
                if all(ismember(PathBNE{k},PathBNE{kk}))
                    PathBNE{k} = [];
                end
            end
        end
        I = cellfun(@(x) ~isempty(x),PathBNE);
        PathBNE = PathBNE(I);
        PathBI = PathBI(I);
        
        if length(PathBI) > 1
            0;
        end
        
        %------------------------------------------------------------------
        % Find point closest to centroid
        %------------------------------------------------------------------
        [I,J] = ind2sub(Size,PathBNE{1});
        mI = mean(I(1:end-1));
        mJ = mean(J(1:end-1));
        dist = sqrt((I - mI).^2 + (J - mJ).^2);
        [~,K] = min(dist);
        
        PathBI1 = PathBI{1}(2:end,1);
        PathJP1 = PathBNE{1};
        JPmin = PathJP1(K);
        PathJP1 = unique(PathJP1);
        PathJP1(PathJP1 == JPmin) = [];
        
        for k = 1 : length(PathJP1)
            BNE(PathBI1(k),:) = [PathJP1(k) JPmin];
        end
        BNE(PathBI1(k+1),:) = [PathJP1(k) JPmin];
    end
end

%--------------------------------------------------------------------------
% Update branch information
%--------------------------------------------------------------------------
BranchNodeEnds(BranchList,:) = BNE;
BranchNodes(BranchList) = mat2cell(BNE,ones(size(BNE,1),1));

%--------------------------------------------------------------------------
% Remove duplicate branches
%--------------------------------------------------------------------------
% [~,I] = unique(BranchNodeEnds,'rows');

n = max(cellfun(@(x) length(x),BranchNodes));
BranchNodes1 = cellfun(@(x) [reshape(x(:),1,length(x)), zeros(1,n - length(x))],BranchNodes,'UniformOutput',0);
BranchNodes1 = vertcat(BranchNodes1{:});
[~,I] = unique(BranchNodes1,'rows');

BranchNodeEnds = BranchNodeEnds(I,:);
BranchNodes = BranchNodes(I);
BranchNodes = cellfun(@(x) x(:),BranchNodes,'UniformOutput',0);


%--------------------------------------------------------------------------
% Merge branches which is splited without joint (not sure why it happens)
%--------------------------------------------------------------------------
BranchNodeEnds = cellfun(@(x) [x(1) x(end)],BranchNodes,'UniformOutput',0);
BranchNodeEnds = cell2mat(BranchNodeEnds(:));
while 1
    for i = 1 : length(BranchNodes)
        [I,J] = find(BranchNodeEnds == BranchNodeEnds(i,1));
        if length(I) == 2
            J = J(~(I == i));
            I = I(~(I == i));
            if isempty(I)
                continue;
            end
            if BranchNodes{i}(1) == BranchNodes{I}(end) %J == 2
                BranchNodes{i} = [BranchNodes{I}(1:end-1); BranchNodes{i}];
                BranchNodes{I} = [];
                
                BranchNodeEnds(i,:) = [BranchNodeEnds(I,1),BranchNodeEnds(i,2)];
                BranchNodeEnds(I,:) = 0;
            elseif BranchNodes{i}(1) == BranchNodes{I}(1)
                BranchNodes{i} = [flip(BranchNodes{I}(2:end)); BranchNodes{i}];
                BranchNodes{I} = [];
                
                BranchNodeEnds(i,:) = [BranchNodeEnds(I,2),BranchNodeEnds(i,2)];
                BranchNodeEnds(I,:) = 0;
            else
                error('Something is wrong');
            end
        end
        
        [I,J] = find(BranchNodeEnds == BranchNodeEnds(i,2));
        if length(I) == 2
            J = J(~(I == i));
            I = I(~(I == i));
            if isempty(I)
                continue;
            end
            if BranchNodes{i}(end) == BranchNodes{I}(1) %J == 1
                BranchNodes{i} = [BranchNodes{i}; BranchNodes{I}(2:end)];
                BranchNodes{I} = [];
                
                BranchNodeEnds(i,:) = [BranchNodeEnds(i,1),BranchNodeEnds(I,2)];
                BranchNodeEnds(I,:) = 0;
            elseif BranchNodes{i}(end) == BranchNodes{I}(end)
                BranchNodes{i} = [BranchNodes{i}; flip(BranchNodes{I}(1:end-1))];
                BranchNodes{I} = [];
                
                BranchNodeEnds(i,:) = [BranchNodeEnds(i,1),BranchNodeEnds(I,1)];
                BranchNodeEnds(I,:) = 0;
            else
                error('Something is wrong');
            end
        end
    end
    id = cellfun(@(x) isempty(x),BranchNodes);
    if nnz(id) == 0
        break;
    end
    BranchNodes(id) = [];
    BranchNodeEnds(id,:) = [];
end

MA_new.Size = Size;
MA_new.BranchNodes = BranchNodes;


