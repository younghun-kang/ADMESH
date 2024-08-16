function id_MA_new = ThinningMA(id_MA,UIFigure)

msg = 'Thinning medial axis...';
uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');

bw = full(id_MA);

%--------------------------------------------------------------------------
% Fill holes with less than 10 pixels
%--------------------------------------------------------------------------
conn = 4;
filled = imfill(bw,conn,'holes');
holes = filled & ~bw;
bigholes = bwareaopen(holes, 10,conn);
smallholes = holes & ~bigholes;
bw = bw | smallholes;
clear filled holes bigholes smallholes;

% %--------------------------------------------------------------------------
% % Expand one pixel (temporary solution to fill holes with small number
% % (around 2~3) of pixels)
% %--------------------------------------------------------------------------
% nhood = [0 1 0; 1 1 1; 0 1 0];
% bw = imdilate(bw,nhood);

%--------------------------------------------------------------------------
% (???) no idea why.. but bwskel looks better than bwmorph...
%--------------------------------------------------------------------------
% id_MA_thinned = bwmorph(bw,'thin');
id_MA_new = bwskel(bw);
% id_MA_thinned = bwmorph(bw,'skel',inf);

%--------------------------------------------------------------------------
% Remove isolated pixels
%--------------------------------------------------------------------------
id_MA_new = bwmorph(id_MA_new,'clean',inf);

id_MA_new = sparse(id_MA_new);
% %--------------------------------------------------------------------------
% % I didn't mean it but this part fills holes
% % (it's better to merge joints latter...)
% % 2021-03-15 This part is commented out as it behaves weird
% %--------------------------------------------------------------------------
% [B_MA_thinned,L_MA_thinned] = bwboundaries(id_MA_thinned);
% nL = zeros(max(L_MA_thinned(:)),1);
% for i = 1 : max(L_MA_thinned(:))
%     nL(i) = nnz(L_MA_thinned == i);
% end
% A_MA_thinned = nL;
% id_MA_thinned = logical(L_MA_thinned > 0);