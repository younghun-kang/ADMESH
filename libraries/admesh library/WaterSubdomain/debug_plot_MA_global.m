
MA_plot = MA_branch;
MA_plot = MA_pruned;
MA_plot = MA_1D;
MA_plot = MA_1D_mainstream;
MA_plot = MA_connected;
% MA_plot = MA;

[I,J] = cellfun(@(x) ind2sub(size(MA_plot.Size),x),MA_plot.BranchNodes,'UniformOutput',0);
temp_x = cellfun(@(x) xg(x),J,'UniformOutput',0);
temp_x = cellfun(@(x) [x(:); nan],temp_x,'UniformOutput',0);
temp_x = vertcat(temp_x{:});
temp_y = cellfun(@(x) yg(x),I,'UniformOutput',0);
temp_y = cellfun(@(x) [x(:); nan],temp_y,'UniformOutput',0);
temp_y = vertcat(temp_y{:});

figure; hold on;
plot(pgon.Vertices([1:end,1],1),pgon.Vertices([1:end,1],2),'k','linewidth',1.5); axis equal;
plot(temp_x,temp_y,'b');

% temp = L2D_boundary;
% temp = cellfun(@(x) [x; nan, nan],temp,'UniformOutput',0);
% temp = vertcat(temp{:});
% plot(temp(:,1),temp(:,2),'b');
% plot(pgon_2D);

%%

figure; hold on;
plot(pgon);
id = AreaInfo.id_MA2D;
[I,J] = ind2sub(size(Dg),id);
plot(xg(J),yg(I),'.');
id = AreaInfo.id_MA1D;
[I,J] = ind2sub(size(Dg),id);
plot(xg(J),yg(I),'.');
axis equal;

%%

figure; hold on;
plot(pgon.Vertices(:,1),pgon.Vertices(:,2),'k'); axis equal;
[I,J] = find(L1D); plot(xg(J),yg(I),'.');
[I,J] = find(L1DtoL2D); plot(xg(J),yg(I),'.');
[I,J] = find(L2D); plot(xg(J),yg(I),'.');










