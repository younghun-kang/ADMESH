function [Mask2DCell,Mask1DCell,M1DtoM2DCell] = FillMask2D(Mask2DCell,Mask1DCell,MA,lfs,delta_filling,delta_A,xg,yg,XY_b,dx,UIFigure)

if ~iscell(Mask2DCell)
    Mask2DCell = {Mask2DCell};
end
if ~iscell(Mask1DCell)
    Mask1DCell = {Mask1DCell};
end
nMask = length(Mask2DCell);

if nMask ~= length(Mask1DCell)
    error;
end

%% ========================================================================
% Identify branch nodes in 2D masks
%==========================================================================
msg = 'Find branch nodes in level 2 masks...';
progdlg = uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg);

nBranch = zeros(nMask,1);
BranchNodes = cell(nMask,1);

for i = 1 : length(MA.BranchNodes)
    id = MA.BranchNodes{i};
    [id1,id2] = ind2sub(size(Mask2DCell{1}),id);
    
    for iMask = 1 : nMask
    I = Mask2DCell{iMask}(id);
    
    if nnz(I) > 0
        I = find(I);
        if any(diff(I) > 1)           
            J = find(diff(I) > 1);
            J = [0; J(:); length(I)];
            J = diff(J);
            I = mat2cell(I(:),J,1);
        else
            I = {I};
        end
        
        for j = 1 : length(I)
            if length(I{j}) > 1
%                 I1 = I{j}([1 end]);
                I1 = I{j}(:);
                I1 = I1(lfs(id(I1)) < delta_filling);
                if ~isempty(I1)
                nBranch(iMask) = nBranch(iMask) + 1;
                BranchNodes{iMask}{nBranch(iMask)} = [id1(I1),id2(I1)];
                end
            end
        end
    end
    end
    
    progdlg.Value = i/length(MA.BranchNodes);
end
% If no branch nodes in some masks, BranchNodes returns empty
IB = find(cellfun(@(x) iscell(x),BranchNodes));
BranchNodes = cellfun(@(x) vertcat(x{:}),BranchNodes(IB),'UniformOutput',0);

%% ========================================================================
% Filling with maximal disk
%==========================================================================
msg = 'Filling with maximal disk...';
progdlg = uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg);

M1DtoM2DCell = cell(nMask,1);
for iMask = 1 : nMask
    if ~ismember(iMask,IB)
        continue;
    end
    x1 = xg(BranchNodes{iMask}(:,2));
    y1 = yg(BranchNodes{iMask}(:,1));
    [~,~,D] = Compute8SSED_v3(XY_b,x1,y1,dx);
    
    [I,J] = find(Mask1DCell{iMask});
    
    K = false(size(I));
    for i = 1 : size(BranchNodes{iMask},1)
        ii = BranchNodes{iMask}(i,1);
        jj = BranchNodes{iMask}(i,2);
        d1 = D(i)/dx;

        I1 = find(abs(I - ii) < d1 & abs(J - jj) < d1);
        d = sqrt((I(I1)-ii).^2 + (J(I1)-jj).^2);
        %     [~,d] = knnsearch([ii,jj],[I,J]);
        K(I1(d < d1)) = 1;
        %     K = sub2ind(size(L1D),I(iii),J(iii));
        progdlg.Value = (iMask-1 + i/size(BranchNodes{iMask},1))/nMask;
    end
    M1DtoM2D = false(size(Mask2DCell{iMask}));
    ii = find(Mask1DCell{iMask});
    M1DtoM2D(ii(K)) = 1;
    
    Mask2DCell{iMask} = (Mask2DCell{iMask} | M1DtoM2D);
    M1DtoM2DCell{iMask} = M1DtoM2D;
end
close(progdlg);


%% ========================================================================
% Remove filled 2D mask regions from 1D masks
%==========================================================================
BranchNodesID = vertcat(MA.BranchNodes{:});
BranchNodesID = unique(BranchNodesID);
for iMask = 1 : nMask
    Mask2D = full(Mask2DCell{iMask});
    Mask1D = full(Mask1DCell{iMask});
    
    Mask1D(M1DtoM2DCell{iMask}) = 0;
    
    %----------------------------------------------------------------------
    % Remove 2D mask not containing any MA points (it may need only when 
    % the filling results isolated mask, which should not happen)
    %----------------------------------------------------------------------
    [I,J] = ind2sub(size(Mask2D),BranchNodesID);
    K = Mask2D(BranchNodesID);
    MA2D = sparse(I,J,K,size(Mask2D,1),size(Mask2D,2));
    
    CC = bwconncomp(full(Mask2D));
    for i = 1 : CC.NumObjects
        id = CC.PixelIdxList{i};
        if ~any(MA2D(id)) || numel(id) < delta_A
            Mask2D(id) = 0;
            Mask1D(id) = 1;
        end
    end
    clear MA2D;
    Mask2DCell{iMask} = Mask2D;
    Mask1DCell{iMask} = Mask1D;
end

%% ========================================================================
% Remove 1D mask regions not containing any MA points 
%==========================================================================    
msg = 'Removing 1D mask regions without MA points...';
progdlg = uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg);
for iMask = 1 : nMask
    Mask2D = full(Mask2DCell{iMask});
    Mask1D = full(Mask1DCell{iMask});
    [I,J] = ind2sub(size(Mask2D),BranchNodesID);
    K = Mask1D(BranchNodesID);
    MA1D = sparse(I,J,K,size(Mask2D,1),size(Mask2D,2));
    
    CC = bwconncomp(full(Mask1D));
    for i = 1 : CC.NumObjects
        id = CC.PixelIdxList{i};
        if ~any(MA1D(id))
            Mask2D(id) = 1;
            Mask1D(id) = 0;
        end
        progdlg.Value = (iMask - 1 + i/CC.NumObjects)/nMask;
    end
    
    Mask2DCell{iMask} = Mask2D;
    Mask1DCell{iMask} = Mask1D;
end

