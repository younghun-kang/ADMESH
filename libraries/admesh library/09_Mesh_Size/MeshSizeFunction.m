function h = MeshSizeFunction(h0,D,hmax,hmin,g,delta,UIFigure)

msg = 'Computing Mesh Size Function...';
uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');

%--------------------------------------------------------------------------
% Initialize Variables
%--------------------------------------------------------------------------
% hn  = hmax*ones(size(h0));
% dt  = delta/2;
% tol = 10e-6;
% R   = 0;
% 
% %--------------------------------------------------------------------------
% % Find all indices that fall within D <= hmax
% %--------------------------------------------------------------------------
% [r,c] = find(D <= hmax/2);
% 
% %--------------------------------------------------------------------------
% % Perform gradient limiting
% %--------------------------------------------------------------------------
% while 1
%     
%     % Compute upwind differences
%     
%     xfordiff    = min((h0(r+1,c) - h0(r,c))./delta,0).^2;
%     
%     xbackdiff   = max((h0(r,c) - h0(r-1,c))./delta,0).^2;
%     
%     yfordiff    = min((h0(r,c+1) - h0(r,c))./delta,0).^2;
%     
%     ybackdiff   = max((h0(r,c) - h0(r,c-1))./delta,0).^2;
%     
%     % Compute Delta
%     Delta = sqrt(xfordiff + xbackdiff + yfordiff + ybackdiff);
%     
%     % Compute next time step
%     hn(r,c) = h0(r,c) + dt*(min(Delta,g) - Delta);
%     
%     % Compute residual
%     R = sum(sum(abs(h0(r,c) - hn(r,c))));
%         
%     if R <= tol; break; end
%     
%     % For next time step
%     h0(r,c) = hn(r,c);
%     
%     drawnow % for graphics
%     
% end

%--------------------------------------------------------------------------
% Perform gradient limiting
%--------------------------------------------------------------------------
h = MeshSizeIterativeSolver(h0,D,hmax,hmin,g,delta);

end