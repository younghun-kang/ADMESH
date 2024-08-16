function MA_1D = ConstructMA1D(MA,Mask1D,Mask2D,UIFigure)

msg = 'Constructing medial axis for open-channel constraints...';
progdlg = uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg);

BranchNodesID = vertcat(MA.BranchNodes{:});
BranchNodesID = unique(BranchNodesID);

I1 = Mask1D(BranchNodesID) | ~Mask2D(BranchNodesID);
nBranch1D = 0;
BranchNodes1D = [];
for i = 1 : length(MA.BranchNodes)
    id = MA.BranchNodes{i};
    
    [~,temp] = ismember(id,BranchNodesID);
    I = I1(temp);
    
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
                nBranch1D = nBranch1D + 1;
                BranchNodes1D{nBranch1D} = id(I{j});
            end
        end
    end
    progdlg.Value = i/length(MA.BranchNodes);
end

MA_1D = [];
MA_1D.Size = MA.Size;
MA_1D.BranchNodes = BranchNodes1D(:);

