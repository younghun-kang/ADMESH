function MA_connected = ConnectMA1Dto2DArea_v2(MA,MaskBoundaryIDs,delta,UIFigure)

msg = 'Connect open-channel/internal boundary constraints to shorelines if they are too close...';
progdlg = uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg);

%% ========================================================================
% Make MA1D is connected to 2D area
%==========================================================================
if isempty(MaskBoundaryIDs) || isempty(MA.BranchNodes)
    MA_connected = MA;
    return;
end

MaskBoundaryIDs_list = vertcat(MaskBoundaryIDs{:});
MaskBoundaryIDs_ind = sub2ind(MA.Size,MaskBoundaryIDs_list(:,1),MaskBoundaryIDs_list(:,2));

BranchNodes = MA.BranchNodes;
BranchNodeEnds = cellfun(@(x) [x(1) x(end)],MA.BranchNodes,'UniformOutput',0);
BranchNodeEnds = cell2mat(BranchNodeEnds(:));

FlagConnected2D = zeros(length(BranchNodes),1);
for i = 1 : length(BranchNodes)
    [x,y] = ind2sub(MA.Size,BranchNodes{i});

    inode = BranchNodes{i};
    
    dist1 = sqrt((MaskBoundaryIDs_list(:,1) - x(1)).^2 + (MaskBoundaryIDs_list(:,2) - y(1)).^2);
    dist2 = sqrt((MaskBoundaryIDs_list(:,1) - x(end)).^2 + (MaskBoundaryIDs_list(:,2) - y(end)).^2);
    
    if nnz(BranchNodeEnds == inode(1)) == 1
    if any(dist1 <= delta & dist1 > 0)
        [~,iend] = min(dist1);
%         MA_connected.ID{i} = vertcat(iend,MA_connected.ID{i});
        x = vertcat(MaskBoundaryIDs_list(iend,1),x);
        y = vertcat(MaskBoundaryIDs_list(iend,2),y);
        inode = vertcat(MaskBoundaryIDs_ind(iend),inode);
        BranchNodeEnds(i,1) = -BranchNodeEnds(i,1);
        FlagConnected2D(i) = FlagConnected2D(i) + 1;
    end
    end
    
    if nnz(BranchNodeEnds == inode(end)) == 1
    if any(dist2 <= delta & dist2 > 0)
        [~,iend] = min(dist2);
%         MA_connected.ID{i} = vertcat(MA_connected.ID{i},iend);
         x = vertcat(x,MaskBoundaryIDs_list(iend,1));
         y = vertcat(y,MaskBoundaryIDs_list(iend,2));
         inode = vertcat(inode,MaskBoundaryIDs_ind(iend));
         BranchNodeEnds(i,2) = -BranchNodeEnds(i,2);
         FlagConnected2D(i) = FlagConnected2D(i) + 1;
    end
    end
%     if nnz(ismember(vertcat(BranchXY{:}),XY([1 end],:),'rows')) > 2
%         FlagConnected2D(i) = FlagConnected2D(i) + 1;
%     end
    
    BranchNodes_connected{i} = inode;
    XY_connected{i} = [x,y];
    Length(i) = sum( sqrt((x(1:end-1) - x(2:end)).^2 + (y(1:end-1) - y(2:end)).^2));
    progdlg.Value = i/length(BranchNodes);
end

n = max(cellfun(@(x) length(x),BranchNodes_connected));
BranchNodes_connected1 = cellfun(@(x) [reshape(x(:),1,length(x)), zeros(1,n - length(x))],BranchNodes_connected,'UniformOutput',0);
BranchNodes_connected1 = vertcat(BranchNodes_connected1{:});
[~,id] = unique(BranchNodes_connected1,'rows');
BranchNodes_connected = BranchNodes_connected(id);
FlagConnected2D = FlagConnected2D(id);

XY_connected = XY_connected(:);
Length = Length(:);

MA_connected.Size = MA.Size;
MA_connected.BranchNodes = BranchNodes_connected;
MA_connected.FlagConnected2D = FlagConnected2D;

