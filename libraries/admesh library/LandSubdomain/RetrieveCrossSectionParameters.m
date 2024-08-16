function CSParams = RetrieveCrossSectionParameters(app,j)


CS = app.CrossSection{j};


m2deg = km2deg(1e-3);
xyzFun = app.xyzFun;
CoordinateSystem = app.CoordinateSystemDropDown.Value;
if strcmpi(CoordinateSystem,'Projected (m)')
    MaxChannelWidth = app.MaxChannelWidthEditField.Value;
    delta = 1;
    unit = 'meter';
elseif strcmpi(CoordinateSystem,'Unprojected (decimal degree)')
    MaxChannelWidth = app.MaxChannelWidthEditField.Value*m2deg;
    delta = 1*m2deg;
    unit = 'degree';
end

cs_length = MaxChannelWidth;
cs_np = ceil(cs_length/delta);
csd = linspace(-cs_length/2,cs_length/2,cs_np)';

CSParams = table;
CSParams.TW(size(CS,1)) = 0;

wbar = waitbar(0);
fwbar = @(x,y) waitbar(x/y,wbar,sprintf('Wairbar (%d/%d)',x,y));
for i = 1 : size(CS,1)
%     if i < 896
%         continue;
%     end
    d = CS.db{i};

    if isempty(d)
        continue;
    end

    csz = CS.fz{i}(csd);
    
    csd1 = csd - CS.cc{i}(1);
    csz1 = csz - CS.cc{i}(2);
    
    csd_center = CS.cc{i}(1);
    
    % Define area and elevation functions of parametric cross-sections
    Acs_p = @(TW,BW,H) (TW+BW)*H/2; % cross-sectional area of parametrized channel
    fTopo_p = @(TW,BW,H,x) (x < -BW/2 & x >= -TW/2).*(-2*H/(TW-BW-1e-6)).*(x + BW/2)...
        + (x > BW/2 & x <= TW/2).*(2*H/(TW-BW-1e-6)).*(x - BW/2);
    
    % Define elevation function of cross-section
    fTopo1 = @(y) interp1(csd1,csz1,y);
    fTopo1_inv = @(y) interp1(csz1,csd1,y);
    
    if nnz(~isnan(csz1)) < 2
        0;
    end

if numel(d) == 1
    twRange = [0 0];
else
    twRange = diff(d)*[.9 1.1];
%     twRange = min(twRange,Width2DElem(i)*.5);
end
%     twRange = [0 .5]*Width2DElem(i);
    
    hRange = [0 0];
    for h = linspace(0,max(csz1),100)
        id1 = find(csz1 > h & csd1 < 0,1,'last');
        id2 = find(csz1 > h & csd1 > 0,1,'first');
        
        if isempty(id1) || isempty(id2)
            continue;
        end
        
        if (csd1(id2) - csd1(id1)) < 0.9*twRange(1)
            hRange(1) = h;
        end
        
        if (csd1(id2) - csd1(id1)) < 1.1*twRange(2)
            hRange(2) = h;
        end
    end
       
    % Start retrieval
    ERR = inf;
    for h = linspace(hRange(1),hRange(2),100)
        id1 = find(csz1 > h & csd1 < 0,1,'last');
        id2 = find(csz1 > h & csd1 > 0,1,'first');
        
        tw = csd1(id2) - csd1(id1);
        
        if isempty(id1) || isempty(id2)
            continue;
        end
        d1 = linspace(csd1(id1),csd1(id2),100);
        fTopo2 = @(y) max(h - fTopo1(y),0);
        Acs = integral(fTopo2,d1(1),d1(end)); % Compute cross-sectional area
        Acs = ComputeCSArea(csd1,csz1,h);
        for bw = linspace(0, tw,100)
            ERR_area = abs(Acs - Acs_p(tw,bw,h))/Acs;
            Topo_p = ParametricChannel(tw,bw,h,d1-mean(d1));
            
            ERR_topo = sqrt(sum( (fTopo1(d1) - Topo_p).^2)/tw);
            
            ERR1 = ERR_area*1 + ERR_topo*1;
            if ERR1 < ERR
                TW = tw;
                BW = bw;
                H = h;
                CSC = mean(d1) + csd_center;
                ERR = ERR1;
            end
        end
    end
    
    if isinf(ERR)
        TW = tw;
        BW = 0;
        H = mean(CS.zb{i} - CS.cc{i}(2));
        CSC = 0;
    end
        
    CSParams.TW (i) = TW;
    CSParams.BW (i) = BW;
    CSParams.H  (i) = H;
    CSParams.CSC(i) = CSC; % cross-section center
    CSParams.ERR(i) = ERR;
    CSParams.BEL(i) = min(csz);
    
    fwbar(i,size(CS,1));
    
end
delete(wbar);

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

id = find(csd2 < csd_center);
[~,i] = min(fxx2(id));

I(1) = id(i);

id = find(csd2 > csd_center);
[~,i] = min(fxx2(id));
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
% figure;
% subplot(211); plot(csd,csz1,csd,ppval(f,csd));
% subplot(212); plot(csd,ppval(fxx,csd));


% idCenterRange = find(abs(csd2) < 2*dx);
% [~,idCenter] = min(ppval(f,csd2(idCenterRange)));
% idCenter = idCenterRange(idCenter);
% csd_center = csd2(idCenter);
csd_center = 0;

N = ceil(max(csz1)/.1);
H = linspace(0.1,max(csz1),N);
K = zeros(size(H));
for j = 1 : length(H)
    h = H(j);
    id1 = find(csz2 > h & csd2 < csd_center,1,'last');
    id2 = find(csz2 > h & csd2 > csd_center,1,'first');

    if ~isempty(id1) && ~isempty(id2)
        d1 = csd2(id1);
        d2 = csd2(id2);
%         d1 = fzero(@(x) fcsz2(x) - h,csd2(id1));
%         d2 = fzero(@(x) fcsz2(x) - h,csd2(id2));
        
        K(j) = sum(ppval(fxx,[d1 d2]));
    else
        K(j) = inf;
    end
end

[~,id] = min(K);
h = H(id);
id1 = find(csz2 > h & csd2 < csd_center,1,'last');
id2 = find(csz2 > h & csd2 > csd_center,1,'first');
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

id = find(csd2 < csd_center);
[~,i] = min(fxx2(id));

I(1) = id(i);

id = find(csd2 > csd_center);
[~,i] = min(fxx2(id));
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




