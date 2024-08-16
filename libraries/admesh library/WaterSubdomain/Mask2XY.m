function [BoundaryXY,BoundaryID] = Mask2XY(Mask,xg,yg,UIFigure)

msg = 'Converting mask boundary to XY coordinates...';
progdlg = uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg);

Mask = full(Mask);
MaskBoundary = bwmorph(Mask,'fill');
MaskBoundary = bwmorph(MaskBoundary,'remove');
MaskBoundary = bwmorph(MaskBoundary,'clean'); % remove isolated points

if nnz(MaskBoundary) == 0
    BoundaryXY = [];
    BoundaryID = [];
    return;
end
k = 0;
BoundaryID = {};
NumBoundaryNodes = nnz(MaskBoundary);
while 1
    [I1,J1] = find(MaskBoundary);
    
    MaskBoundary1 = [];
    FSTEP = {'N','E','S','W'};
    DIR = {'clockwise','clockwise','counterclockwise','counterclockwise'};
    for i = 1 : 4
        try 
            MaskBoundary1 = bwtraceboundary(Mask,[I1(1),J1(1)],FSTEP{i},8,inf,DIR{i});
            break;
        catch
            msg = sprintf('Boundary search failed in ''%s'' direction.\n',DIR{i});
            uiconfirm(UIFigure,msg,'ADMESH',...
                'Options',{'OK'},'DefaultOption',1,'Icon','Error');
            
        end
    end

    if isempty(MaskBoundary1)
        error('Boundary serach failed in all directions.');
    end
    
    %----------------------------------------------------------------------
    % Add boundary segment if no duplicated segments
    %----------------------------------------------------------------------
    if ~any(cellfun(@(x) all(ismember(x,MaskBoundary1,'rows')),BoundaryID))
        k = k + 1;
        BoundaryID{k} = MaskBoundary1;
%         W2D_boundary{k} = sub2ind(size(W2D_bw_boundary),w2dbdry(:,1),w2dbdry(:,2));
    end
        
    %----------------------------------------------------------------------
    % Remove indexes of added boundary segment
    %----------------------------------------------------------------------
    IdRemove = sub2ind(size(MaskBoundary),MaskBoundary1(:,1),MaskBoundary1(:,2));
    MaskBoundary(IdRemove) = 0;
    
    if nnz(MaskBoundary) == 0
        break;
    end
    
%     IdRemove = find(ismember([I1,J1],w2dbdry,'rows'));
%     I1(IdRemove) = [];
%     J1(IdRemove) = [];
    
%     if isempty(I1)
%         break;
%     end
    progdlg.Value = (1 - length(I1)/NumBoundaryNodes);
end
clear W2D_bw_boundary;

% figure; hold on;
% for i = 1 : length(W2D_boundary)
%     plot(W2D_boundary{i}(:,2),W2D_boundary{i}(:,1));
% end

%--------------------------------------------------------------------------
% Remove boundary with only one vertex
%--------------------------------------------------------------------------
BoundaryID = BoundaryID(:);
I = (cellfun(@(x) size(unique(x,'rows'),1) == 1,BoundaryID));
BoundaryID(I) = [];

%--------------------------------------------------------------------------
% Convert indexes to xy-coordinates
%--------------------------------------------------------------------------
BoundaryXY = cellfun(@(x) [xg(x(:,2)), yg(x(:,1))],BoundaryID,'UniformOutput',false);
BoundaryXY = BoundaryXY(:);
