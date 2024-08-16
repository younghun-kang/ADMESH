function [Vx,Vy] = Compute8SSED(Xq,Yq,Bx,By,dx)

Debugging = 0;
%--------------------------------------------------------------------------
% Interpolate boundary points so that the resolution matchs to the
% background grid
%--------------------------------------------------------------------------
t = 1 : length(Bx);
t2 = 1;
dist = sqrt((Bx(1:end-1)-Bx(2:end)).^2 + (By(1:end-1)-By(2:end)).^2);
if isempty(dx)
    dx = mean(mean(diff(Xq,1,2)));
end
for i = 1 : length(Bx)-1
    n = max(2,round(dist(i)/dx));
    t2 = [t2, linspace(t2(end),t2(end)+1,n)];
end
t2 = unique(t2);
t2 = t2(:);

Bx2 = interp1(t,Bx,t2);
By2 = interp1(t,By,t2);

id2boundary = knnsearch([Bx2,By2],[Xq(:),Yq(:)],'k',1);

Vx = Bx2(id2boundary(:,1)) - Xq(:);
Vy = By2(id2boundary(:,1)) - Yq(:);
Vx = reshape(Vx,size(Xq,1),size(Xq,2));
Vy = reshape(Vy,size(Xq,1),size(Xq,2));

if Debugging == 1
    figure; hold on;
    plot(Bx,By,'k');
    quiver(Xq,Yq,Vx,Vy);
end