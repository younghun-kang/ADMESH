function iSK = Skeletonization_DivSSED_v2(BX,BY,Vx,Vy)

Dibugging = 0;

DIV = divergence(BX,BY,Vx,Vy);

W = sqrt(Vx.^2 + Vy.^2);
Vx1 = Vx./W;
Vy1 = Vy./W;
DIV1 = divergence(BX,BY,Vx1,Vy1);

iSK = DIV > 0;
% iSK = DIV./W > 1e4;

if Dibugging == 1
% id = 1 : 10 : numel(BX);
% id = id(:);
% id = iSK;
% figure; hold on;
% quiver(BX(id),BY(id),U1(id),V1(id));
% axis equal tight;

figure; hold on;
% surf(BX,BY,DIV,'edgecolor','none'); view(2);
myScatter3(BX(iSK),BY(iSK),W(iSK));
axis equal tight;
title('||V||');

figure; hold on;
% surf(BX,BY,DIV,'edgecolor','none'); view(2);
myScatter3(BX(iSK),BY(iSK),DIV(iSK));
axis equal tight;
title('Div');

figure; hold on;
% surf(BX,BY,DIV,'edgecolor','none'); view(2);
myScatter3(BX(iSK),BY(iSK),DIV(iSK)./W(iSK));
axis equal tight;
title('DIV/||V||');

figure; hold on;
surf(BX,BY,DIV,'edgecolor','none');
axis equal tight;
view(2);
end




