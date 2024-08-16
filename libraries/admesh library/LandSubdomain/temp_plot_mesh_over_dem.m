figure;
AXf = axes;
AXf(2) = axes;

axes(AXf(2));
myPlot14File(MBRW_fine,'field','mesh');
myLinkAxes(AXf);

figure;
AXc = axes;
AXc(2) = axes;

axes(AXc(2));
myPlot14File(MBRW_coarse,'field','mesh');
myLinkAxes(AXc);

%%
ax = axis;

[I,J] = find(X10 > ax(1) & X10 < ax(2) & Y10 > ax(3) & Y10 < ax(4));
I = unique(I);
J = unique(J);

axes(AXf(1));
surf(X10(I,J),Y10(I,J),Z10(I,J),'EdgeColor','none','FaceColor','interp');
view(2); axis equal tight; mycmap('dem'); colorbar;
axis(ax);
axes(AXf(2));


axes(AXc(1));
surf(X10(I,J),Y10(I,J),Z10(I,J),'EdgeColor','none','FaceColor','interp');
view(2); axis equal tight; mycmap('dem'); colorbar;
axis(ax);
axes(AXc(2));


%%
figure;
AX_dem = axes;
AX_dem(2) = axes;

axes(AX_dem(2));
% myPlot14File(MBRW,'field','mesh');
hold on;
% plot(PARAMS.Txr,PARAMS.Tyr,'k');
% plot(PARAMS.Txl,PARAMS.Tyl,'k');

plot(PARAMS_bank.Bxr,PARAMS_bank.Byr,'k');
plot(PARAMS_bank.Bxl,PARAMS_bank.Byl,'k');
% plot(PARAMS_banksym.Bxr,PARAMS_banksym.Byr,'k');
% plot(PARAMS_banksym.Bxl,PARAMS_banksym.Byl,'k');
%%
% figure(1);
% AX_dem = AX_sym;
ax = axis;
cax = caxis(AX_dem(1));
[I,J] = find(X > ax(1) & X < ax(2) & Y > ax(3) & Y < ax(4));
I = unique(I);
J = unique(J);

surf(AX_dem(1),X(I,J),Y(I,J),Z(I,J),'EdgeColor','none','FaceColor','interp');
view(AX_dem(1),2);
axis(AX_dem(1),'equal');
colorbar(AX_dem(1));
if isequal(cax,[0 1])
    cax = zlim(AX_dem(1));
end
mycmap('dem',AX_dem(1),cax);

myLinkAxes(AX_dem);
% 

% SyncFigs([5 4],'cax');











