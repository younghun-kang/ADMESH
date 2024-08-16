% 
% tic; 
% for i = 1 : 10
%     [D,gradD] = SignedDistanceFunction_v2(PTS,X,Y,delta,h_max);
% end
% toc;
% tic; 
% for i = 1 : 10
% [Vx,Vy] = Compute8SSED(X,Y,x,y,delta/2);
% end
% toc;
% tic; 
% for i = 1 : 10
% [Vx,Vy] = Compute8SSED_v2(x,y,X,Y);
% end
% toc;
% 
% tic; 
% for i = 1 : 10
% [Vx,Vy] = Compute8SSED_v3(XYb,X,Y,delta);
% end
% toc;


[D,gradD] = SignedDistanceFunction_v2(PTS,X,Y,delta,hmax);
[Vx,Vy] = Compute8SSED_v3(XY_b,X,Y,delta);