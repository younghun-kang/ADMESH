function [Mask2D,Mask1D] = SplitMask(Mask,lfs,delta_lfs,delta_A,MA)

ID = vertcat(MA.BranchNodes{:});
ID = unique(ID);
[I,J] = ind2sub(size(Mask),ID);
MaskMA = sparse(I,J,true(size(I)),size(Mask,1),size(Mask,2));

Mask2D = full(Mask & abs(lfs) >= delta_lfs+1e-8);
MA2D = full(MaskMA & Mask2D);

CC = bwconncomp(Mask2D);
% Mask2D = false(size(Mask));
for i = 1 : CC.NumObjects
    id = CC.PixelIdxList{i};
    if numel(id) < delta_A || ~any(MA2D(id))
        Mask2D(id) = 0;
    end
end
% L2D = sparse(L2D);
clear MASK MA2D;

Mask1D = full(Mask);
Mask1D(Mask2D) = 0;

%% ========================================================================
%==========================================================================
% Mask1D = full(Mask & abs(lfs) < delta_lfs);
% Mask2D1 = Mask;
% Mask2D1(Mask1D) = 0;
% MA2D = full(MaskMA & Mask2D1);
% 
% CC = bwconncomp(Mask2D1);
% Mask2D = false(size(Mask));
% for i = 1 : CC.NumObjects
%     id = CC.PixelIdxList{i};
%     if any(MA2D(id)) && numel(id) > delta_A
%         Mask2D(id) = 1;
%     end
% end
% % L2D = sparse(L2D);
% clear MASK MA2D;
