function MA_new = MergeShortBridge(MA,MinChannelLength,dx)

Size = MA.Size;
BranchNodes = MA.BranchNodes;

BranchNodeEnds = cellfun(@(x) [x(1) x(end)],BranchNodes,'UniformOutput',0);
BranchNodeEnds = cell2mat(BranchNodeEnds(:));

if isfield(MA,'BranchXY')
    BranchXY = MA.BranchXY;
end
if isfield(MA,'FlagConnected2D')
    FlagConnected2D = MA.FlagConnected2D;
end
while 1
    [I,J] = cellfun(@(x) ind2sub(Size,x),BranchNodes,'UniformOutput',0);
    BranchLength = cellfun(@(x,y) sum(sqrt(diff(x).^2 + diff(y).^2))*dx,I,J);
    
    for i = 1 : size(BranchNodeEnds,1)
        if BranchLength(i) > MinChannelLength
            continue;
        end
        

        [I,J] = ind2sub(Size,BranchNodes{i});
        
        mI = mean(I);
        mJ = mean(J);
        dist = sqrt((I - mI).^2 + (J - mJ).^2);
        [~,K] = min(dist);
        K = sub2ind(Size,I(K),J(K));
        
        [I1,J1] = find(BranchNodeEnds == BranchNodeEnds(i,1));
        J1 = J1(I1 ~= i);
        I1 = I1(I1 ~= i);
        
        [I2,J2] = find(BranchNodeEnds == BranchNodeEnds(i,2));
        J2 = J2(I2 ~= i);
        I2 = I2(I2 ~= i);
        
        if length(I1) < 2 || length(I2) < 2
            continue;
        end
        
        for k = 1 : length(I1)
            if J1(k) == 1
                BranchNodes{I1(k)} = [K; BranchNodes{I1(k)}(2:end)];
            else 
                BranchNodes{I1(k)} = [BranchNodes{I1(k)}(1:end-1); K];
            end 
        end 
        
        for k = 1 : length(I2)
            if J2(k) == 1
                BranchNodes{I2(k)} = [K; BranchNodes{I2(k)}(2:end)];
            else 
                BranchNodes{I2(k)} = [BranchNodes{I2(k)}(1:end-1); K];
            end 
        end 
        
        BranchNodes{i} = [];
        BranchNodeEnds(i,:) = 0;
        
        
%         if length(I) == 2
%             J = J(~(I == i));
%             I = I(~(I == i));
%             
%             if J == 1
%                 JointConnectivity(i,:) = [JointConnectivity(i,1),JointConnectivity(I,2)];
%                 BranchNodes{i} = [BranchNodes{i}; BranchNodes{I}];
%                 BranchNodes{I} = [];
%             else
%                 JointConnectivity(i,:) = [JointConnectivity(i,1),JointConnectivity(I,1)];
%                 BranchNodes{i} = [BranchNodes{i}; flip(BranchNodes{I})];
%                 BranchNodes{I} = [];
%             end
%         end
    end
    id = cellfun(@(x) isempty(x),BranchNodes);
    if nnz(id) == 0
        break;
    end
    BranchNodes(id) = [];
    BranchNodeEnds(id,:) = [];
    if isfield(MA,'BranchXY')
        BranchXY(id) = [];
    end
    if isfield(MA,'FlagConnected2D')
    FlagConnected2D(id) = [];
    end
end

MA_new.Size = Size;
MA_new.BranchNodes = BranchNodes;
