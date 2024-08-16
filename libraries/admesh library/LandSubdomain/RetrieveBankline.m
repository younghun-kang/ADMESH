function Params = RetrieveBankline(fDEM,LineXY,MaxChannelWidth,cs_delta,Method,varargin)
%==========================================================================
% Parser optional arguments
%==========================================================================
%--------------------------------------------------------------------------
% Set default values
%--------------------------------------------------------------------------
Unit = 'meter';
%--------------------------------------------------------------------------
% Get optional values
%--------------------------------------------------------------------------
SkipNext = false;
for i = 1 : length(varargin)
    if SkipNext
        SkipNext = false;
        continue;
    end
    ivar = varargin{i};
    if any(strcmpi(ivar,{'meter','degree'}))
        Unit = ivar;
    else
        error('Something is wrong with optional arguments');
    end
end

%==========================================================================
Params = table;
Params.d{size(LineXY,1)} = 0;
if numel(MaxChannelWidth) == 1
    MaxChannelWidth = repmat(MaxChannelWidth,size(LineXY,1),1);
end

for i = 1 : size(LineXY,1)
%     if i < 1617
%         continue;
%     end
    XYq = LineXY(i,:);
    
    j = knnsearch(LineXY,XYq);
    
    px = LineXY(j,1);
    py = LineXY(j,2);
    if i == 1
        tx = LineXY(j+1,1) - LineXY(j,1);
        ty = LineXY(j+1,2) - LineXY(j,2);
    elseif i == size(LineXY,1)
        tx = LineXY(j,1) - LineXY(j-1,1);
        ty = LineXY(j,2) - LineXY(j-1,2);
    else
        tx = LineXY(j+1,1) - LineXY(j-1,1);
        ty = LineXY(j+1,2) - LineXY(j-1,2);
    end
    nx = ty/sqrt(tx^2 + ty^2);
    ny = -tx/sqrt(tx^2 + ty^2);
    
%     twRange = TWRange(i,:);
    cs_length = 1*MaxChannelWidth(i);
    cs_np = ceil(cs_length/cs_delta);
    csd = linspace(-cs_length/2,cs_length/2,cs_np)';
    csx = px + nx*csd;
    csy = py + ny*csd;
    
    % Linear interpolation.
    csz = fDEM(csx,csy);
    csz = double(csz);
    
    if strcmpi(Unit,'degree')
        csd1 = deg2km(csd)*1e3;
    elseif strcmpi(Unit,'meter')
        csd1 = csd;
    end
    
    I = find(islocalmin(csz));
    if ~isempty(I)
        csd_center = csd(I);
        [~,isort] = sort(abs(csd_center));
        csd_center = csd_center(isort(1));
        csz_min = csz(I(isort(1)));
    else
        csd_center = 0;
        csz_min = fDEM(px,py);        
    end

    csd1 = csd1 - csd_center;

    % Set the bottom of cross-sectional elevation is zero
    csz1 = csz - csz_min;
    if nnz(~isnan(csz1)) < 2
        0;
    end
    
    csd1 = csd1(~isnan(csz1));
    csz1 = csz1(~isnan(csz1));
    
    csd1 = double(csd1);
    csz1 = double(csz1);
    % Set channel top-width based on curvature

    SmoothingRMSE = .01*max(csz1);
    [csz0,SmoothingParam] = csapsRMSE(csd1,csz1,SmoothingRMSE,csd1);
    csz0 = csz0 + csz_min;
    
%     SmoothingParam = 1e-2;
    
    if Method == 1
    [d,z] = CurvatureBasedEstimation(csd1,csz1,SmoothingParam);
    elseif Method == 2
        [d,z] = CurvatureBasedEstimation2(csd1,csz1,SmoothingParam);
    elseif Method == 2.1
        [d,z] = CurvatureBasedEstimation21(csd1,csz1,SmoothingParam);
    elseif Method == 3
        [d,z] = CurvatureBasedEstimation3(csd1,csz1,SmoothingParam);
    elseif Method == 4
        [d,z] = CurvatureBasedEstimation4(csd1,csz1,SmoothingParam);
    elseif Method == 5
        [d,z] = CurvatureBasedEstimation5(csd1,csz1,SmoothingParam);
    elseif Method == 6
        [d,z] = CurvatureBasedEstimation6(csd1,csz1);
    end
    
    d = d + csd_center;
    if strcmpi(Unit,'degree')
        d = km2deg(d*1e-3);
        csd1 = km2deg(csd1*1e-3);
        csd_center = km2deg(csd_center*1e-3);
    end
    
    Params.d{i} = d;
    Params.z{i} = z+csz_min;
    Params.c{i} = [csd_center,csz_min];
   
end

end

function z = ParametricChannel(TW,BW,H,x)

if TW > BW
    z = (x < -BW/2 & x >= -TW/2).*(-2*H/(TW-BW)).*(x + BW/2)...
        + (x > BW/2 & x <= TW/2).*(2*H/(TW-BW)).*(x - BW/2)...
        + (abs(x) > TW/2).*H + (abs(x) <= BW/2)*.0;
else
    z = (abs(x) > TW/2).*H + (abs(x) <= BW/2)*.0;
end

end

function Acs = ComputeCSArea(csd1,csz1,h)
fTopo1 = @(d) interp1(csd1,csz1,d);
fTopo2 = @(d) max(h - fTopo1(d),0);

id1 = find(csz1 > h & csd1 < 0,1,'last');
id2 = find(csz1 > h & csd1 > 0,1,'first');
if isempty(id1)
    id1 = 1;
end
if isempty(id2)
    id2 = length(csd1);
end

if ~isempty(id1) && ~isempty(id2)
    d1 = linspace(csd1(id1),csd1(id2),100);
    Acs = integral(fTopo2,d1(1),d1(end)); % Compute cross-sectional area
else
    Acs = nan;
end
end

function [d,z] = CurvatureBasedEstimation(csd1,csz1,SmoothingParam)
% dx: resolution of DEM 
%% ========================================================================
% Retrieve cross section parameters using curvature
%==========================================================================
csd2 = linspace(csd1(1),csd1(end),2*length(csd1));
csz2 = interp1(csd1,csz1,csd2);

% f = interp1(csd,double(csz1),'spline','pp');
f = csaps(csd1,double(csz1),SmoothingParam);
fx = fnder(f,1);
fxx = fnder(f,2);

csd_center = 0;

f2 = ppval(f,csd2);
fx2 = ppval(fx,csd2);
fxx2 = ppval(fxx,csd2);

F = @(x) ppval(fxx,x)./(1+ppval(fx,x).^2).^(3/2);

id = find(csd2 < csd_center);
[~,i] = min(F(id));

I(1) = id(i);

id = find(csd2 > csd_center);
[~,i] = min(F(id));
I(2) = id(i);

d = csd2(I);
z = csz2(I);


%     figure;
%     plot(csd,csz1,'--','color',[.5 .5 .5]); hold on;
%     plot(csd2,ppval(f,csd2),'b');
%     % set(gca,'PlotBoxAspectRatio',[2.5 1 1]);
%     hold on; plot(csd2(icurvmin),ppval(f,csd2(icurvmin)),'rs');
%     hold on; plot(csd2(icurvmax),ppval(f,csd2(icurvmax)),'bs');
%     hold on; plot(csd2(idmin),ppval(f,csd2(idmin)),'ks','MarkerFaceColor','k');
%     hold on; plot(csd2([ibankl ibankr]),ppval(f,csd2([ibankl ibankr])),'rs','markerfacecolor','r');
%     hold on; plot(csd2([icurveminl2,icurveminr2]),f2([icurveminl2,icurveminr2]),'vb','markerfacecolor','b');
%     % hold on; plot(csd(id),fo(csd(id)),'ok','MarkerfaceColor','k');
%     grid on;
%     % axis(ax);
%     xlabel('Distance (m)'); ylabel('Topography (m)');
%======================================================================

end

function [d,z] = CurvatureBasedEstimation2(csd1,csz1,SmoothingParam)
%% ========================================================================
% Retrieve cross section parameters using curvature
%==========================================================================
csd2 = linspace(csd1(1),csd1(end),2*length(csd1));
csz2 = interp1(csd1,csz1,csd2);
fcsz2 = @(d) interp1(csd1,csz1,d);
% f = interp1(csd,double(csz1),'spline','pp');
f = csaps(csd1,double(csz1),SmoothingParam);
fx = fnder(f,1);
fxx = fnder(f,2);

f2 = ppval(f,csd2);
fx2 = ppval(fx,csd2);
fxx2 = ppval(fxx,csd2);

% [~,id] = min(K);
% h = H(id);
% id1 = find(csz2 > h & csd2 < csd_center,1,'last');
% id2 = find(csz2 > h & csd2 > csd_center,1,'first');

F = @(x) ppval(fxx,x)./(1+ppval(fx,x).^2).^(3/2);
% F = @(x) ppval(fxx,x) + abs(ppval(fx,x));
% F = @(x) ppval(fx,x)/max(abs(fx2))+ppval(fxx,x)/max(abs(fxx2));

I = find(csd2 < 0);
[~,id] = max(ppval(f,csd2(I)));
d_maxL = csd2(I(id));
I = find(csd2 > 0);
[~,id] = max(ppval(f,csd2(I)));
d_maxR = csd2(I(id));


minH = 0.1;
I = find(csd2 < -5 & f2 > minH);
if isempty(I)
    I = find(csd2 < 0 & f2 > minH);
    if isempty(I)
        I = find(csd2 < 0);
    end
end
[~,id1] = min(F(csd2(I)));
id1 = I(id1);

I = find(csd2 > 5 & f2 > minH);
if isempty(I)
    I = find(csd2 > 0 & f2 > minH);
    if isempty(I)
        I = find(csd2 > 0);
    end
end
[~,id2] = min(F(csd2(I)));
id2 = I(id2);

d1 = csd2(id1);
d2 = csd2(id2);
% d1 = fzero(@(x) fcsz2(x) - h,csd2(id1));
% d2 = fzero(@(x) fcsz2(x) - h,csd2(id2));
d = [d1 d2];
z = fcsz2([d1 d2]);

% figure; plot(csd,csz1);
% hold on; plot(d,z,'o');

end

function [d,z] = CurvatureBasedEstimation21(csd1,csz1,SmoothingParam)
%% ========================================================================
% Retrieve cross section parameters using curvature
%==========================================================================
csd2 = linspace(csd1(1),csd1(end),2*length(csd1));
csz2 = interp1(csd1,csz1,csd2);
fcsz2 = @(d) interp1(csd1,csz1,d);
% f = interp1(csd,double(csz1),'spline','pp');
f = csaps(csd1,double(csz1),SmoothingParam);
fx = fnder(f,1);
fxx = fnder(f,2);

f2 = ppval(f,csd2);
fx2 = ppval(fx,csd2);
fxx2 = ppval(fxx,csd2);
% figure;
% subplot(211); plot(csd,csz1,csd,ppval(f,csd));
% subplot(212); plot(csd,ppval(fxx,csd));

                  
            
                                                    
                                                     

% idCenterRange = find(abs(csd2) < 2*dx);
% [~,idCenter] = min(ppval(f,csd2(idCenterRange)));
% idCenter = idCenterRange(idCenter);
% csd_center = csd2(idCenter);
csd_center = 0;

dH = 0.01;
minH = 0.1;
N = ceil(max(csz1)/dH);
H = linspace(0,max(csz1),N);
F = @(x) ppval(fxx,x)./(1+ppval(fx,x).^2).^(3/2);
K = zeros(size(H));
for j = 1 : length(H)
    h = H(j);
    id1 = find(csz2 >= h & csd2 < csd_center,1,'last');
    id2 = find(csz2 >= h & csd2 > csd_center,1,'first');
                          
             

                               
    if ~isempty(id1) && ~isempty(id2)
        d1 = csd2(id1);
        d2 = csd2(id2);
%         d1 = fzero(@(x) fcsz2(x) - h,csd2(id1));
%         d2 = fzero(@(x) fcsz2(x) - h,csd2(id2));
        
        K(j) = sum(F([d1 d2]));
    else
        K(j) = inf;
    end
end
                          
             

[~,id] = min(K);
h = H(id);
id1 = find(csz2 >= h & csd2 < csd_center,1,'last');
id2 = find(csz2 >= h & csd2 > csd_center,1,'first');
if isempty(id1) || isempty(id2)
    0;
end
d1 = csd2(id1);
d2 = csd2(id2);
% d1 = fzero(@(x) fcsz2(x) - h,csd2(id1));
% d2 = fzero(@(x) fcsz2(x) - h,csd2(id2));
d = [d1 d2];
z = fcsz2([d1 d2]);

% figure; plot(csd,csz1);
% hold on; plot(d,z,'o');

end



function [d,z] = CurvatureBasedEstimation3(csd1,csz1,SmoothingParam)
% dx: resolution of DEM 
%% ========================================================================
% Retrieve cross section parameters using curvature
%==========================================================================
csd2 = linspace(csd1(1),csd1(end),2*length(csd1));
csz2 = interp1(csd1,csz1,csd2);

% f = interp1(csd,double(csz1),'spline','pp');
f = csaps(csd1,double(csz1),SmoothingParam);
fx = fnder(f,1);
fxx = fnder(f,2);

csd_center = 0;

f2 = ppval(f,csd2);
fx2 = ppval(fx,csd2);
fxx2 = ppval(fxx,csd2);
F = @(x) ppval(fxx,x)./(1+ppval(fx,x).^2).^(3/2);
K = F(csd2);

id = find(csd2 < csd_center);
[~,i] = min(K(id));

I(1) = id(i);

id = find(csd2 > csd_center);
[~,i] = min(K(id));
I(2) = id(i);

d = csd2(I);
z = csz2(I);

h = min(csz2(I));
id1 = find(csz1 >= h & csd1 < csd_center,1,'last');
id2 = find(csz1 >= h & csd1 > csd_center,1,'first');
if isempty(id1) || isempty(id2)
    0;
end
d = csd1([id1 id2]);
z = csz1([id1 id2]);

end

function [d,z] = CurvatureBasedEstimation4(csd1,csz1,SmoothingParam)
%% ========================================================================
% Retrieve cross section parameters using curvature
%==========================================================================
csd2 = linspace(csd1(1),csd1(end),2*length(csd1));
csz2 = interp1(csd1,csz1,csd2);
fcsz2 = @(d) interp1(csd1,csz1,d);
% f = interp1(csd,double(csz1),'spline','pp');
f = csaps(csd1,double(csz1),SmoothingParam);
fx = fnder(f,1);
fxx = fnder(f,2);

f2 = ppval(f,csd2);
fx2 = ppval(fx,csd2);
fxx2 = ppval(fxx,csd2);
% figure;
% subplot(211); plot(csd,csz1,csd,ppval(f,csd));
% subplot(212); plot(csd,ppval(fxx,csd));


% idCenterRange = find(abs(csd2) < 2*dx);
% [~,idCenter] = min(ppval(f,csd2(idCenterRange)));
% idCenter = idCenterRange(idCenter);
% csd_center = csd2(idCenter);
csd_center = 0;

dH = 0.01;
minH = 0.1;
N = ceil(max(csz1)/dH);
H = linspace(0,max(csz1),N);
id_left = zeros(size(H));
id_right = zeros(size(H));
for j = 1 : length(H)
    h = H(j);
    id = find(csz2 >= h & csd2 < csd_center,1,'last');
    if ~isempty(id)
        id_left(j) = id;
    end
    
    id = find(csz2 >= h & csd2 > csd_center,1,'first');
    if ~isempty(id)
        id_right(j) = id;
    end
end

I = id_left > 0 & id_right > 0;
id_left = id_left(I);
id_right = id_right(I);

d_left = csd2(id_left);
d_right = csd2(id_right);

H1 = H(I);
F = @(x) ppval(fxx,x)./(1+ppval(fx,x).^2).^(3/2);
K = F(d_left) + F(d_right);

D = d_right - d_left;

I = find(H1 > minH & D > 10);
if isempty(I)
    I = find(H1 > minH);
    if isempty(I)
        I = 1:length(K);
    end
end
[~,id] = min(K(I));
d = [d_left(I(id)), d_right(I(id))];
z = fcsz2(d);

end

function [d,z] = CurvatureBasedEstimation5(csd1,csz1,SmoothingParam)
%% ========================================================================
% Retrieve cross section parameters using curvature
%==========================================================================
csd2 = linspace(csd1(1),csd1(end),2*length(csd1));
csz2 = interp1(csd1,csz1,csd2);
fcsz2 = @(d) interp1(csd1,csz1,d);
% f = interp1(csd,double(csz1),'spline','pp');
f = csaps(csd1,double(csz1),SmoothingParam);
fx = fnder(f,1);
fxx = fnder(f,2);

f2 = ppval(f,csd2);
fx2 = ppval(fx,csd2);
fxx2 = ppval(fxx,csd2);
% figure;
% subplot(211); plot(csd,csz1,csd,ppval(f,csd));
% subplot(212); plot(csd,ppval(fxx,csd));


% idCenterRange = find(abs(csd2) < 2*dx);
% [~,idCenter] = min(ppval(f,csd2(idCenterRange)));
% idCenter = idCenterRange(idCenter);
% csd_center = csd2(idCenter);
csd_center = 0;

dH = 0.01;
minH = 0.1;
N = ceil(max(csz1)/dH);
H = linspace(0,max(csz1),N);
H = H(:);
id_left = zeros(size(H));
id_right = zeros(size(H));
for j = 1 : length(H)
    h = H(j);
    id = find(csz2 >= h & csd2 < csd_center,1,'last');
    if ~isempty(id)
        id_left(j) = id;
    end
    
    id = find(csz2 >= h & csd2 > csd_center,1,'first');
    if ~isempty(id)
        id_right(j) = id;
    end
end

idvalid_left = (id_left > 0);
idvalid_right = (id_right > 0);

d_left = csd2(id_left(idvalid_left));
d_right = csd2(id_right(idvalid_right));

H_left = H(idvalid_left);
H_right = H(idvalid_right);

F = @(x) ppval(fxx,x)./(1+ppval(fx,x).^2).^(3/2);
K_left = F(d_left);
K_right = F(d_right);

K_left = K_left(:)';
K_right = K_right(:);
K = K_left + K_right;

d_left = d_left(:)';
d_right = d_right(:);
D = d_right - d_left;

kk = 1 : length(K(:));
kk = kk(:);
[ii,jj] = ind2sub(size(K),kk);

I = true(size(ii));
I1 = I & abs(H_right(ii)-H_left(jj))./((H_right(ii)+H_left(jj))/2) < 0.1;
if nnz(I1) > 0
    I = I1;
end
I1 = I & H_right(ii) > minH & H_left(jj) > minH;
if nnz(I1) > 0
    I = I1;
end
I1 = I & D(:) > 10;
if nnz(I1) > 0
    I = I1;
end


kk = sub2ind(size(K),ii(I),jj(I));

[~,id] = min(K(kk));

[ii,jj] = ind2sub(size(K),kk(id));

d = [d_left(jj), d_right(ii)];
z = [H_left(jj), H_right(ii)];


end             


function [d,z] = CurvatureBasedEstimation6(csd1,csz1)
%% ========================================================================
% Retrieve cross section parameters using curvature
%==========================================================================
csd2 = linspace(csd1(1),csd1(end),1*length(csd1));
csz2 = interp1(csd1,csz1,csd2);

csd_center = 0;

dH = 0.01;
minH = 0.1;
N = ceil(max(csz1)/dH);
H = linspace(0,max(csz1),N);
H = H(:);
id_left = zeros(size(H));
id_right = zeros(size(H));
for j = 1 : length(H)
    h = H(j);
    id = find(csz2 >= h & csd2 < csd_center,1,'last');
    if ~isempty(id)
        id_left(j) = id;
    end
    
    id = find(csz2 >= h & csd2 > csd_center,1,'first');
    if ~isempty(id)
        id_right(j) = id;
    end
end

W = zeros(size(H));
W2 = zeros(size(H));
for j = 1 : length(H)
    if id_left(j) > 0 && id_right(j) > 0
        W2(j) = csd2(id_right(j)) - csd2(id_left(j));
    end
end

dW = diff(W)./diff(H);
dW2 = diff(W2)./diff(H);
% d = 0;
% z = 0;
idvalid = find(id_left > 0 & id_right > 0);
idvalid(idvalid == length(W)) = [];
idvalid(H(idvalid) < minH) = [];
[~,id] = max(dW2(idvalid));
id = idvalid(id);

d = [csd2(id_left(id)), csd2(id_right(id))];
z = [H(id), H(id)];


end    

