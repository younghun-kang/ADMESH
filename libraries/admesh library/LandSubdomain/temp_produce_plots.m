
prefix = 'figs/comp_single_vs_bankline/';
postfix = '_zoom1-3';

myexport_fig(figure(1),[prefix,'smooth_single_param',postfix]);
myexport_fig(figure(2),[prefix,'smooth_bankline_mesh',postfix]);
myexport_fig(figure(3),[prefix,'smooth_single_dem',postfix]);
myexport_fig(figure(4),[prefix,'smooth_bankline_dem',postfix]);


%%

AX_temp = {AX_ss_dem,AX_sb_dem};
ax = axis(AX_temp{1}(1));
cax = caxis(AX_temp{1}(1));


[I,J] = find(X > ax(1) & X < ax(2) & Y > ax(3) & Y < ax(4));
I = unique(I);
J = unique(J);

for i = 1 : length(AX_temp)
    AX = AX_temp{i};
    surf(AX(1),X(I,J),Y(I,J),Z(I,J),'EdgeColor','none','FaceColor','interp');
    view(AX(1),2);
    axis(AX(1),ax);
    colorbar(AX(1));
    mycmap('dem',AX(1),cax);
end
 
AX_temp = {AX_ss_dem,AX_ss_param,AX_sb_dem,AX_sb_mesh};
for i = 1 : length(AX_temp)
    AX = AX_temp{i};
    view(AX(1),2);
    axis(AX(1),ax);
    colorbar(AX(1));
    mycmap('dem',AX(1),cax);
end

colorbar(AX_sb_mesh,'off');

myLinkAxes({AX_ss_dem,AX_ss_param,AX_sb_dem,AX_sb_mesh},'Visible','Clim');
