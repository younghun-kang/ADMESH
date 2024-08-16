
MA_plot = MA_branch;
MA_plot = MA_pruned;
MA_plot = MA_1D;
MA_plot = MA_1D_mainstream;
MA_plot = MA_connected;
% MA_plot = MA;

figure; hold on;
plot(x_b,y_b,'k');
for i = 1 : MA_plot.nBranch
    id = MA_plot.BranchNodes{i};
    plot(X(id),Y(id),'linewidth',1.2);
%     plot3(MA_plot.XY{i}(:,1),MA_plot.XY{i}(:,2),i*ones(size(MA_plot.XY{i},1),1),'linewidth',1.2);

%     myPlot3(X(id),Y(id),Vmag(id) + D(id),'linewidth',1.2);
end
axis equal;


%%
% nhood = [0 1 0; 1 1 1; 0 1 0];
% bw = imdilate(id_MA_original,nhood);
figure; hold on; axis equal; myc = mycolors;
plot(x_b,y_b,'k');
% plot(X(bw),Y(bw),'o','color',myc.or);
plot(X(id_MA_original),Y(id_MA_original),'.r');
% plot(X(id_MA_thinned),Y(id_MA_thinned),'.b');
% plot(X(id_MA1D_Filled),Y(id_MA1D_Filled),'.b');


%%
Z = AreaInfo.L2D;
Z = AreaInfo_Filled.L2D;
figure; hold on;
surf(X,Y,Z,'edgeColor','none');
myPlotOnTop(x_b,y_b,'k');
axis equal; mycmap('jet',max(Z(:))); colorbar;