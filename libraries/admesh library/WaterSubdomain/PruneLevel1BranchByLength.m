function MA = PruneLevel1BranchByLength(dx,MA,MinChannelLength,FlagPreserve)

BranchNodeEnd = cellfun(@(x) [x(1) x(end)],MA.BranchNodes,'UniformOutput',0);
BranchNodeEnd = cell2mat(BranchNodeEnd(:));
BranchNodeEnd = abs(BranchNodeEnd);
BranchNodes = MA.BranchNodes;
if isfield(MA,'BranchXY')
    BranchXY = MA.BranchXY;
end
if ~exist('FlagPreserve','var')
    FlagPreserve = false(length(BranchNodes),1);
else
%     FlagPreserve = FlagPreserve >= 2;
end
% MinChannelLength = 2/(DEG2KM*1e3);

while 1
    
    %----------------------------------------------------------------------
    % Pruning short branch
    %----------------------------------------------------------------------
    [I,J] = cellfun(@(x) ind2sub(MA.Size,x),BranchNodes,'UniformOutput',0);
    BranchLength = cellfun(@(x,y) sum(sqrt(diff(x).^2 + diff(y).^2))*dx,I,J);
    id = find(BranchLength < MinChannelLength);
    
    for i = 1 : length(id)
        j = id(i);
        num_nghb(1) = nnz(BranchNodeEnd == BranchNodeEnd(j,1));
        num_nghb(2) = nnz(BranchNodeEnd == BranchNodeEnd(j,2));
        if ~any(num_nghb == 1)
            id(i) = 0;
        end
        if max(num_nghb) > 1 && FlagPreserve(j) == 1
            id(i) = 0;
        end
    end
    
    id = nonzeros(id);
    id(FlagPreserve(id) == 2) = [];
    
    if isempty(id)
        break;
    end
    
    BranchNodes(id) = [];
    BranchNodeEnd(id,:) = [];
    if isfield(MA,'BranchXY')
        BranchXY(id) = [];
    end
    FlagPreserve(id) = [];
end
I = cellfun(@(x) ~isempty(x),BranchNodes);
BranchNodes = BranchNodes(I);
FlagPreserve = FlagPreserve(I);

n = max(cellfun(@(x) length(x),BranchNodes));
BranchNodes1 = cellfun(@(x) [reshape(x(:),1,length(x)), zeros(1,n - length(x))],BranchNodes,'UniformOutput',0);
BranchNodes1 = vertcat(BranchNodes1{:});
[~,I] = unique(BranchNodes1,'rows');

BranchNodes = BranchNodes(I);
FlagPreserve = FlagPreserve(I);

MA.BranchNodes = BranchNodes;
MA.FlagConnected2D = FlagPreserve;



