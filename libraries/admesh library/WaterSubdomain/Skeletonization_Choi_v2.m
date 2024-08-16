function iSK = Skeletonization_Choi_v2(BX,BY,Vx,Vy,rho,rho2,dtheta)

Debugging = 0;

N = size(BX,1);
M = size(BX,2);

% rho = 1e-2;
iSK = zeros(N,M);
dsq = zeros(numel(iSK),1);
theta1 = zeros(numel(iSK),1);
for j = 2 : M-1
    for i = 2 : N-1
        k = N*(j-1) + i;
%         if ~ismember(k,idIN)
%             continue;
%         end
%         k_nghb = k + [-M-1, -M, -M+1, -1, 1, M-1, M, M+1];
        k_nghb = k + [-N-1, -N, -N+1, -1, 1, N-1, N, N+1];
        for ii = 1 : 8
            kk = k_nghb(ii);
            
            Dsq = (Vx(kk) - Vx(k))^2 + (Vy(kk) - Vy(k))^2;
            COSINE = (Vx(kk)*Vx(k) + Vy(kk)*Vy(k))...
                /(sqrt(Vx(kk)^2 + Vy(kk)^2) * sqrt(Vx(k)^2 + Vy(k)^2));
            if isnan(COSINE)
                COSINE = 1;
            end
            
            COSINE = min(COSINE,1);
            COSINE = max(COSINE,-1);
            
            theta = acos(COSINE);
            
            Qi_norm = Vx(kk)^2 + Vy(kk)^2;
            Q_norm = Vx(k)^2 + Vy(k)^2;
            
            if Dsq >= rho^2 && (Qi_norm - Q_norm) <= max(BX(k) - BX(kk),BY(k)-BY(kk))
                dsq(k) = Dsq; theta1(k) = theta;
                iSK(k) = 1;
            end
        end
        %------------------------------------------------------------------
        % rejection corners
        %------------------------------------------------------------------
%         if sqrt(dsq(k)) < rho2  && theta1(k) <= dtheta
%             iSK(k) = 0;
%         end
    end
end
% id_sk = reshape(id_sk,N,M);
iSK = logical(iSK);

if Debugging == 1
figure;
hold on; myScatter3(BX(iSK),BY(iSK),dsq(iSK));
caxis([0 rho^2]);
title('Dsq');

% figure;
% hold on; myScatter3(BX(iSK),BY(iSK),-D(iSK));
% caxis([0 rho]);
% title('D');

figure;
hold on; myScatter3(BX(iSK),BY(iSK),theta1(iSK));
% quiver(BX(iSK),BY(iSK),SSED(iSK,1),SSED(iSK,2),10,'color','k');
caxis([0 dtheta]);
title('angle');

id = 1 : 2 : numel(BX);
id = id(:);
figure; hold on;
quiver(BX(id),BY(id),SSED(id,1),SSED(id,2),10);
axis equal tight;

end




