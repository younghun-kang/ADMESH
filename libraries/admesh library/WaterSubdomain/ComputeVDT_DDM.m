function [Vx,Vy,Mask1,Mask2] = ComputeVDT_DDM(pgon,dd_ID,dd_Boxes,xg,yg,delta,UIFigure)

D = sparse(length(yg),length(xg));
Vx = sparse(length(yg),length(xg));
Vy = sparse(length(yg),length(xg));
Mask1 = false(length(yg),length(xg));
Mask2 = false(length(yg),length(xg));

msg = 'Computing Vector Distance Transform...';
progdlg = uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');

for iDD = 1 : size(dd_ID,1)
       
    x = xg(dd_ID{iDD,2});
    y = yg(dd_ID{iDD,1});
    [X,Y] = meshgrid(x,y);
    
    Dx = max(x) - min(x);
    Dy = max(y) - min(y);
    
    x1 = min(x) - Dx;
    x2 = max(x) + Dx;
    y1 = min(y) - Dy;
    y2 = max(y) + Dy;
    
    nghbbox = polyshape([x1 x2 x2 x1],[y1 y1 y2 y2]);
    
    ipgon_nghb = intersect(pgon,nghbbox);
    ipgon_nghb_regions = regions(ipgon_nghb);
    A = area(ipgon_nghb_regions);
    [~,I] = sort(A,'descend');
    ipgon_nghb_regions = ipgon_nghb_regions(I);
    
    
    k1 = 0;
    XY_b1 = [];
    for iRR = 1 : length(ipgon_nghb_regions)
        iregion = ipgon_nghb_regions(iRR);
        % iregion = ipgon_nghb;
        
        temp = intersect(iregion,dd_Boxes(iDD));
        if temp.NumRegions == 0
            continue;
        end
        iregion_noholes = rmholes(iregion);
        I = find(isnan(iregion.Vertices(:,1)));
        I = [0; I(:); size(iregion.Vertices,1)+1];
        k1 = k1 + 1;
        XY_b1{k1} = iregion_noholes.Vertices([1:end,1],:); % boundary including the dd-box
        
        for i = 1 : length(I)-1
            i1 = I(i) + 1;
            i2 = I(i+1) - 1;
            J = [i1:i2, i1];
            
            if all(cellfun(@(x) nnz(~ismember(x,iregion.Vertices(J,:))),XY_b1) > 0)
                k1 = k1 + 1;
                XY_b1{k1} = iregion.Vertices(J,:);
            end
        end
    end
    XY_b1{end+1} = nghbbox.Vertices([1:end,1],:);
    
    %--------------------------------------------------------------------------
    % Compute distance from shoreline and gradient of it
    %--------------------------------------------------------------------------
    PTS = [];
    PTS.Poly(1).x = XY_b1{1}(:,1);
    PTS.Poly(1).y = XY_b1{1}(:,2);
    for i = 2 : length(XY_b1)
        PTS.Poly(i).x = XY_b1{i}(:,1);
        PTS.Poly(i).y = XY_b1{i}(:,2);
    end
    % [D,gradD] = SignedDistanceFunction_v2(PTS,X,Y,delta,hmax);
    % % D = round(D,16);
    
    IN = PointsInDomain3(X,Y,PTS);
    IN = ~IN; % Flip since it is closed by land boundary (nghbbox)
    if nnz(IN) == 0
        continue;
    end
    
    %--------------------------------------------------------------------------
    % Compute 8SSED map
    %--------------------------------------------------------------------------
    % [Vx,Vy] = Compute8SSED(X,Y,x_b,y_b,delta);
    [Vxloc,Vyloc] = Compute8SSED_v3(XY_b1,X,Y,delta);
    % [Vx,Vy,D] = Compute8SSED_v4(ipgon,X,Y,delta);
%     Vx(~IN) = 0; Vy(~IN) = 0;
    Vmloc = sqrt(Vxloc.^2 + Vyloc.^2);
    
%     Dloc(IN) = -Dloc(IN);
%     Dg1 = full(D(dd_ID{iDD,1},dd_ID{iDD,2}));
%     I = abs(Dg1) == 0 | abs(Dg1) > abs(Dloc); % negative values for outside domain
% %     I = D < Dg1; % ignore outside domain
%     Dg1(I) = Dloc(I);
    
%     D(dd_ID{iDD,1},dd_ID{iDD,2}) = Dg1;
    Vx1 = full(Vx(dd_ID{iDD,1},dd_ID{iDD,2}));
    Vy1 = full(Vy(dd_ID{iDD,1},dd_ID{iDD,2}));
    Vm1 = sqrt( Vx1.^2 + Vy1.^2);
    I = Vm1 == 0 | ( Vm1 > Vmloc);
    
    Vx1(I) = Vxloc(I);
    Vy1(I) = Vyloc(I);
    Vx(dd_ID{iDD,1},dd_ID{iDD,2}) = Vx1;
    Vy(dd_ID{iDD,1},dd_ID{iDD,2}) = Vy1;
    
    
    Mask1(dd_ID{iDD,1},dd_ID{iDD,2}) = IN;
    Mask2(dd_ID{iDD,1},dd_ID{iDD,2}) = ~IN;
    
%     clear D Dg1 Vx Vy Vxg1 Vyg1 Vmag Vmag1 X Y;
    progdlg.Indeterminate = 'off';
    progdlg.Value = iDD/size(dd_ID,1);
end

