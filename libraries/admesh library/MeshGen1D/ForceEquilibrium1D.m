function [p,P,ConvCriteria] = ForceEquilibrium1D(xStart,xEnd,fh,h0,fixedPoints)

%--------------------------------------------------------------------------
% Set constant parameters
%--------------------------------------------------------------------------
Fscale = 1.2; DT = .2;
MaxIter = 1e5;
ConvTol = 1e-4;

%--------------------------------------------------------------------------
% Create (uniform) initial distribution
%--------------------------------------------------------------------------
pInit = xStart : h0 : xEnd;
pInit = pInit(:);

%--------------------------------------------------------------------------
% Apply rejection method and add fixed points
%--------------------------------------------------------------------------
r0 = 1./(fh(pInit).^2);
pInit = pInit(rand(length(pInit),1) < r0/max(r0));
pInit = [fixedPoints(:); pInit];
pInit = unique(pInit);
pInit = sort(pInit);
if isequal(pInit,fixedPoints)
    p = pInit;
    P{1} = p;
    ConvCriteria = 0;
    return;
end
%--------------------------------------------------------------------------
% Preallocating arrays for storing results
%--------------------------------------------------------------------------
P = cell(MaxIter,1); 
ConvCriteria = zeros(MaxIter,1);

%--------------------------------------------------------------------------
% Start force equilibrium
%--------------------------------------------------------------------------
p = pInit;
P{1} = p;
nIter = 0;
% wbar = waitbar(0);
% fwbar = @(x,y) waitbar(x/y,wbar,sprintf('Force equilibrium 1D... (%.2f)',x/y));
while 1
%     nIter = nIter + 1;
    pMid = (p(2:end) + p(1:end-1))/2;
    pMid = pMid(:);
    
    %----------------------------------------------------------------------
    % Compute h-values between points by taking mean value (of 11 points)
    %----------------------------------------------------------------------
    W1 = 0:.1:1;
    W2 = 1:-.1:0;
    
    p1 = p(1:end-1);
    p2 = p(2:end);
%     hMid = mean(fh(p1*W1 + p2*W2),2);
%     hMid = max(h(p1*W1 + p2*W2),[],2);
    hMid = fh((p1+p2)/2);
%     hMid = .5*mean(h(p1*W1 + p2*W2),2) + .5*h((p1+p2)/2);

    %----------------------------------------------------------------------
    % Compute element lengths
    %----------------------------------------------------------------------
    L_now = p(2:end) - p(1:end-1);
    
    L_desired = hMid.*Fscale.*sqrt(sum(L_now.^2)/sum(hMid.^2));
%     L_desired = hMid;
%     L_desired(L_desired < h_min) = h_min;
%     L_desired(L_desired > h_max) = h_max;
    
    %----------------------------------------------------------------------
    % Compute forces on each element
    %----------------------------------------------------------------------
    F = max(L_desired - L_now,0);
%     F = L_desired - L_now;
    Ftot = [0; F(1:end-1) - F(2:end); 0];
    
    %----------------------------------------------------------------------
    % Zero out forces on fixed points
    %----------------------------------------------------------------------
    idfixed = (ismember(p,fixedPoints));
    Ftot(idfixed) = 0;
%     if all(Ftot == 0)
%         n10 = round(length(p)/10);
%         p = [p(:); xEnd*rand(n10,1)];
%         continue;
%     end

    %----------------------------------------------------------------------
    % Update points locations
    %----------------------------------------------------------------------
    p = p + DT*Ftot;
    p = sort(p);
    
    %----------------------------------------------------------------------
    % Bring back outside points
    %----------------------------------------------------------------------
    if any(p < xStart)
        id = p < xStart;
        p(id) = xStart + fh(xStart);
    end
    if any(p > xEnd)
        id = p > xEnd;
        p(id) = xEnd - fh(xEnd);
        %     break;
    end
    p = sort(p);

    %----------------------------------------------------------------------
    % Save updated points
    %----------------------------------------------------------------------
    nIter = nIter + 1;
    P{nIter+1} = p;
    
    %----------------------------------------------------------------------
    % Check tolerence and maximum iteration number
    %----------------------------------------------------------------------
    ConvCriteria(nIter) = max(abs(DT*Ftot(2:end-1))/h0);
    if ConvCriteria(nIter) < ConvTol
        break;
    end
    
    if nIter > MaxIter
        [~,id] = min(ConvCriteria);
        p = P{id};
        warning([...
            'Potential issues are in ForceEquilibrium1D program. ',...
            'Program ended as maximum iteration reached. ',...
            'Mesh points with the least movement is chosen as output.']);
        break;
    end
%      if ConvTol/ConvCriteria(nIter) > nIter/MaxIter
%          fwbar(ConvTol,ConvCriteria(nIter));
%      else
%          fwbar(nIter,MaxIter);
%      end
    
end
% delete(wbar);

P(nIter+2 : MaxIter) = [];
ConvCriteria(nIter+1 : MaxIter) = [];


end