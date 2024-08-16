function L2D_new_boundary = GetLayerBoundary_v2(L,dx)
%==========================================================================
% GetLayerBoundary_v2
% - This function does this, this, and this...
% 
% Update history
% ????-??-?? (v1) - written by Younghun Kang
% 2021-03-15 (v2) - output all boundaries with 'Smoothing' option if bwboundaries identifies more than one boundary
% 
%==========================================================================

% L : layer
% X,Y background grid
% BX,BY : boundary points
%

Debugging = 0;
nL = max(L(:));
L2D_new_boundary = cell(nL,1);
dy = dx;
dx = max(dx,dy);

iB = 0;
for iL = 1 : nL
    iL2D = (L == iL);
    
        B = bwboundaries(double(iL2D));
        if length(B) > 1
            warning(...
                ['Multiple boundaries found for layer #%d. '...
                'Proceed with assumping that they correspond holes. '...
                'Otherwise, it might result wrong.'],iL);
        end
        
        bx = []; by = [];
        for i = 1 : length(B)
            x = B{i}(:,2);
            y = B{i}(:,1);
            
%             %--------------------------------------------------------------
%             % Approach 1: fix start and end point of boundary (it does not
%             % guarantee continuity of derivatives)
%             %--------------------------------------------------------------            
%             SmoothingParam = 0.1;
%             id_fixedPoints = [1,length(x)];
%             t = 1 : length(x);
%             w = ones(length(x),1);
%             w(id_fixedPoints) = 1e8;
%             Smooth_x = csaps(t,x,SmoothingParam,t,w);
%             Smooth_y = csaps(t,y,SmoothingParam,t,w);
%             Smooth_x(id_fixedPoints) = x(id_fixedPoints);
%             Smooth_y(id_fixedPoints) = y(id_fixedPoints);
            
            %--------------------------------------------------------------
            % Approach 2: wind-up boundary twice and get middle portion of
            % it, i.e., extent [0,1] to [0,2] and get [1/2,3/2]
            %--------------------------------------------------------------
            x = [x; x];
            y = [y; y];
            SmoothingParam = 0.1;
            id_fixedPoints = [];
            t = 1 : length(x);
            w = ones(length(x),1);
            w(id_fixedPoints) = 1e8;
            Smooth_x = csaps(t,x,SmoothingParam,t,w);
            Smooth_y = csaps(t,y,SmoothingParam,t,w);
            
            n = round(length(x)/4);
            Smooth_x = Smooth_x(n:3*n);
            Smooth_y = Smooth_y(n:3*n);

            
            bx = Smooth_x(:);
            by = Smooth_y(:);
            
%             bx = ppval(csape(1:length(bx),bx,'periodic'),linspace(1,length(bx),length(bx)*2));
%             by = ppval(csape(1:length(by),by,'periodic'),linspace(1,length(by),length(by)*2));
            
            bxy = unique([bx(:),by(:)],'rows','stable');
            if size(bxy,1) < 2
                continue;
            end
            bxy(end+1,:) = bxy(1,:);
            
            iB = iB + 1;
            L2D_new_boundary{iB} = bxy;
        end
    
%     L2D_new_boundary{iL} = [bx(:),by(:)];
end


if Debugging == 1
    figure; hold on; myc = mycolors; daspect([1 1 1]);
%     plot(x_b,y_b,'k');
%     myPlotOnTop(x,y,'.k');
    plot(Smooth_x,Smooth_y,'r');
    plot(bx,by,'b');
    plot(Smooth_x(end),Smooth_y(end),'ob');
    plot(bxy(:,1),bxy(:,2),'color',myc.gn);
%     plot(L2D_new_boundary(:,1),L2D_new_boundary(:,2),'b')
    axis equal; box on;
end














