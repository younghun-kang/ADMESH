function [dd_pgons,dd_Boxes,dd_ID,xg,yg] = DecomposePolyshape(pgon,Nx,Ny,dx,dBox,UIFigure)

msg = 'Decomposing input mask...';
progdlg = uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');

minx = min(pgon.Vertices(:,1));
maxx = max(pgon.Vertices(:,1));
miny = min(pgon.Vertices(:,2));
maxy = max(pgon.Vertices(:,2));

DDwidth = dBox;
DDheight = dBox;

DDx = linspace(minx - DDwidth/2,maxx + DDwidth/2,Nx+1);
DDy = linspace(miny - DDheight/2,maxy + DDheight/2,Ny+1);

xg = min(DDx)-dx : dx : max(DDx)+dx;
yg = min(DDy)-dx : dx : max(DDy)+dx;
xg = xg(:);
yg = yg(:);

k = 0;
dd_pgons = polyshape;
dd_Boxes = polyshape;
for i = 1 : Nx
    for j = 1 : Ny
        boxx = DDx([i i+1 i+1 i]);
        boxx([1 4]) = boxx([1 4]) - dx;
        boxx([2 3]) = boxx([2 3]) + dx;
        
        boxy = DDy([j j j+1 j+1]);
        boxy([1 2]) = boxy([1 2]) - dx;
        boxy([3 4]) = boxy([3 4]) + dx;
        
        ibox = polyshape(boxx, boxy);

        ipgon = intersect(pgon,ibox);
        
        if ipgon.NumRegions > 0
            k = k + 1;
            dd_pgons(k) = ipgon;
            dd_Boxes(k) = ibox;
        end
        n1 = Nx*(i-1)+j;
        n2 = Nx*Ny;
        progdlg.Indeterminate = 'off';
        progdlg.Value = n1/n2;
    end
end



dd_ID = cell(length(dd_pgons),2);
for i = 1 : length(dd_pgons)
    ipgon = dd_pgons(i);
    
    I1 = find(yg <= min(ipgon.Vertices(:,2) - DDheight/2),1,'last');
    I2 = find(yg >= max(ipgon.Vertices(:,2) + DDheight/2),1,'first');
    II = I1:I2;
    
    J1 = find(xg <= min(ipgon.Vertices(:,1) - DDwidth/2),1,'last');
    J2 = find(xg >= max(ipgon.Vertices(:,1)) + DDwidth/2,1,'first');
    JJ = J1:J2;
    
    dd_ID{i,1} = II(:);
    dd_ID{i,2} = JJ(:);
end