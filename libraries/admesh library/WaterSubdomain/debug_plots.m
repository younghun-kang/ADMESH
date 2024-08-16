
figure; hold on;
plot(x_b,y_b,'k');
% plot(X(id_MA),Y(id_MA),'.r');
% plot(X(bw),Y(bw),'.r');
plot(X(id_MA_thinned),Y(id_MA_thinned),'.b');

axis equal;

temp = bwmorph(id_MA_thinned,'branchpoints');
plot(X(temp),Y(temp),'ob');
% plot(X(MA_branch.ID_Joint),Y(MA_branch.ID_Joint),'or');

return;
%%
bw = id_MA;
[B_MA_thinned,L_MA_thinned] = bwboundaries(bw);
figure; hold on;
plot(x_b,y_b,'k');
% plot(X(id_MA),Y(id_MA),'.r');
for i = 1 : max(L_MA_thinned(:))
    plot(X(L_MA_thinned==i),Y(L_MA_thinned==i),'.');
end
% plot(X(id_MA_thinned),Y(id_MA_thinned),'.b');

axis equal;



bw = id_MA;
id_hole = find(bw == 1);
N = size(bw,1); M = size(bw,2);
id_hole(id_hole <= N | id_hole > (M-1)*N) = [];
bw(id_hole+1) = 1;
bw(id_hole-1) = 1;
bw(id_hole+N) = 1;
bw(id_hole-N) = 1;

figure; hold on; axis equal;
plot(x_b,y_b,'k');
% plot(X(id_MA),Y(id_MA),'.r');
plot(X(bw),Y(bw),'.');

