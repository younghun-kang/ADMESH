function L2D_new_boundary = GetLayerBoundary(L,X,Y,x_b,y_b,option)
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
dx = mean(mean(diff(X,1,2)));
dy = mean(mean(diff(Y,1,1)));
dx = max(dx,dy);

iB = 0;
for iL = 1 : nL
    iL2D = (L == iL);
    if strcmpi(option,'Fill-Cut')
        t = 1 : length(x_b);
        t2 = 1;
        dist = sqrt((x_b(1:end-1)-x_b(2:end)).^2 + (y_b(1:end-1)-y_b(2:end)).^2);
        
        for i = 1 : length(x_b)-1
            n = max(2,round(dist(i)/dx));
            t2 = [t2, linspace(t2(end),t2(end)+1,n)];
        end
        t2 = unique(t2);
        t2 = t2(:);
        
        x_b2 = interp1(t,x_b,t2);
        y_b2 = interp1(t,y_b,t2);
        
        x = X(iL2D); y = Y(iL2D);
        [~,dist] = knnsearch([x,y],[x_b2,y_b2]);
        I = dist < sqrt(2)*dx;
        x1 = vertcat(X(iL2D),x_b2((I)));
        y1 = vertcat(Y(iL2D),y_b2((I)));
        iB = boundary(x1,y1,1);
        bx = x1(iB); by = y1(iB);
        
        iB = iB + 1;
        L2D_new_boundary{iB} = [bx(:),by(:)];
    elseif strcmpi(option,'Smoothing')
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
            
            % Project boundary points to current coordinate system
            b = size(X,2); a = 1;
            d = max(X(:)); c = min(X(:));
            x = (d-c)/(b-a)*(x - a) + c;
            b = size(Y,1); a = 1;
            d = max(Y(:)); c = min(Y(:));
            y = (d-c)/(b-a)*(y - a) + c;
            
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














