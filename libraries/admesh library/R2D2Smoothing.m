function R2D2Smoothing(app,nConnectingLines)

msg = 'Applying R2D2 Smoothing...';
uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');

app.xyzFun_old = app.xyzFun;

Constraints = app.PTS.Constraints;
if isempty(Constraints)
    msg = 'No internal constraint found.';
    uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','warning');
    return;
end
I = abs([Constraints.num]) == 18;
FL = {Constraints(I).xy};

FL2 = cellfun(@(x) [x; nan(1,2)],FL,'UniformOutput',0);
FL2 = vertcat(FL2{:});

% Set end points of each channels
basinXY = cellfun(@(x) [x(end,:)],FL,'UniformOutput',0);
basinXY = vertcat(basinXY{:});

% Add second end point of channels if the end point is a joint point
for i = 1 : size(basinXY,1)
    if nnz(ismember(FL2,basinXY(i,:),'rows')) > 1
        basinXY(end+1,:) = FL{i}(end-1,:);
    end
end

% % try add all points
% basinXY = cellfun(@(x) [x(:,:)],FL,'UniformOutput',0);
% basinXY = vertcat(basinXY{:});

% Find draniagebasins using TopoToolbox
FD = app.TTC.FD;
L = drainagebasins(FD,basinXY(:,1),basinXY(:,2));

% Set BasinID array (need check)
BasinID = double(fliplr(L.Z'));
BasinID(BasinID == 0) = NaN;

% Check number of basins
nBasins = max(BasinID(:));
if nBasins ~= size(basinXY,1)
    error('Error in checking number of basins.');
end

% Check number of id for each basin
nID = zeros(nBasins,1);
for i = 1 : length(nID)
    nID(i) = nnz(BasinID(:) == i);
end

xyzFun = app.xyzFun;

x = xyzFun.GridVectors{1};
y = xyzFun.GridVectors{2};
[X,Y] = ndgrid(x,y);
Z = xyzFun.Values;

% Visualize
% figure; hold on;
% x = xyzFun.GridVectors{1};
% y = xyzFun.GridVectors{2};
% [X,Y] = ndgrid(x,y);
% Z = xyzFun.Values;
% 
% surf(X,Y,BasinID,'edgecolor','none'); view(2);
% % plot3(BP(:,1),BP(:,2),1e5*ones(size(BP(:,2))),'k','linewidth',2);
% plot3(FL2(:,1),FL2(:,2),nBasins*ones(size(FL2(:,1))),'b','linewidth',2);
% axis equal tight;
% colormap('turbo');
% colorbar;
% 
% label_x = zeros(nBasins,1);
% label_y = zeros(nBasins,1);
% label_text = split(num2str(1:nBasins));
% for i = 1 : nBasins
%     label_x(i) = mean(X(BasinID == i));
%     label_y(i) = mean(Y(BasinID == i));
% end
% 
% text(label_x,label_y,nBasins*ones(size(label_x)),label_text);

%%

Znew = NaN(size(Z));

% nConnectingLines = 10;

tic;
for iBasin = 1 : nBasins
    %% Find subbasin boundary and corresponding channel
    K = find(BasinID == iBasin);
    [I,J] = ind2sub(size(BasinID),K);
    I = min(I):max(I);
    J = min(J):max(J);

    BasinMask = false(size(BasinID));
    BasinMask(BasinID == iBasin) = 1;

    % Find subbasin boundary
    B = bwboundaries(BasinID == iBasin);
    % Remove redundant boundary (not sure why it happens)
    I = cellfun(@(x) size(x,1) > 2,B);
    B = B(I);
    % check boundary
    if length(B) ~= 1
        error('Multiple subbasin boundaries are extracted.');
        figure; hold on;
        for j = 1 : length(B)
            plot(B{j}(:,1),B{j}(:,2));
        end
    end
    B = B{1};

    % Find subbasin boundary coordinates
    Bx = x(B(:,1));
    By = y(B(:,2));
    Bx = Bx(:);
    By = By(:);

    % Find channel points within subbasin 
    temp = zeros(length(FL),1);
    for j = 1 : length(FL)
        x1 = FL{j}(:,1);
        y1 = FL{j}(:,2);
        [inS,onS] = inpolygon(x1,y1,Bx,By);
        temp(j) = nnz(inS) + nnz(onS);
    end

    % Find channel id corresponding to the subbasin
    [~,j] = max(temp);

    % Find channel coordinates
    x1 = FL{j}(:,1);
    y1 = FL{j}(:,2);

    % Check if any point of channel is out of subbasin
    [inS,onS] = inpolygon(x1,y1,Bx,By);
    I = find(inS | onS);
    if max(diff(I)) > 1
        warning('Some of channel points are out of subbasin');
    end

    % Set channel points
    Cx = x1(I);
    Cy = y1(I);
    Cz = -xyzFun(Cx,Cy);
    
    % Find subbasin boundary point closest to the end of channel
    iStart = knnsearch([Bx,By],[Cx(end),Cy(end)]);
    
    % Reorder subbasin boundary
    if iStart ~= 1
        Bx = Bx([iStart:end,1:iStart-1]);
        By = By([iStart:end,1:iStart-1]);
    end
    Bz = -xyzFun(Bx,By);

    %% Method 3: use triangulation
    %----------------------------------------------------------------------
    % Set points for triangulation
    %----------------------------------------------------------------------
    P1 = [Bx, By];
    P2 = [Cx, Cy];

    I = ismember(P2,P1,'rows');
    P2(I,:) = [];
    
    P1 = unique(P1,'rows','stable');
    P2 = unique(P2,'rows','stable');
    P = [P1; P2];
    %----------------------------------------------------------------------
    % Check duplicated points
    %----------------------------------------------------------------------
%     P = unique(P,'rows','stable');
% 
%     I = ismember(P1,P,'rows');
%     P1 = P1(I,:);
% 
%     I = ismember(P2,P,'rows');
%     P2 = P2(I,:);

    %----------------------------------------------------------------------
    % Construct triangulation
    %----------------------------------------------------------------------
    TRI = delaunayTriangulation(P);
    x1 = P(:,1);
    y1 = P(:,2);

%     figure; triplot(TRI);

    %----------------------------------------------------------------------
    % Retrieve connectivity list
    %----------------------------------------------------------------------
    T = TRI.ConnectivityList;

%     % Reset points
%     P = TRI.Points;
%     
%     I = ismember(P,P1,'rows');
%     P1 = P(I,:);
% 
%     I = ismember(P,P2,'rows');
%     P2 = P(I,:);

    %----------------------------------------------------------------------
    % Remove triangulation out of subbasin (half side)
    %----------------------------------------------------------------------
    mx = mean(x1(T),2);
    my = mean(y1(T),2);

    [inS,onS] = inpolygon(mx,my,P1(:,1),P1(:,2));
    I = ~inS & ~onS;
    T(I,:) = [];

%     figure; triplot(T,P(:,1),P(:,2),'k'); axis equal tight;

    %----------------------------------------------------------------------
    % Find elements containing one boundary node and one channel node
    %----------------------------------------------------------------------
    xT = x1(T);
    yT = y1(T);

    xT = xT';
    yT = yT';

    % Check if edge contains boundary vertices 
    I1 = ismember([xT(:), yT(:)],P1,'rows');
    I1 = reshape(I1,size(xT,1),size(xT,2));

    % Check if edge contains channel vertices
    I2 = ismember([xT(:), yT(:)],P2,'rows');
    I2 = reshape(I2,size(xT,1),size(xT,2));

    % Find element containing only boundary or channel vertices 
    I = sum(I1).*sum(I2) == 0;
    T(I,:) = [];
    I1(:,I) = [];
    I2(:,I) = [];

%     hold on; triplot(T,P(:,1),P(:,2),'b'); axis equal tight;

    %----------------------------------------------------------------------
    % Find edges connecting channel and boundary
    %----------------------------------------------------------------------
    % Find edge from boundary to channels
    I3 = I1 & I2([2 3 1],:);

    [II1,J] = find(I3);
    II2 = II1+1;
    II2(II2 == 4) = 1;

    K1 = sub2ind(size(T),J,II1);
    K2 = sub2ind(size(T),J,II2);
    C1 = [T(K1),T(K2)];
    
    % Find edge from channel to boundary
    I3 = I1 & I2([3 1 2],:);

    [II1,J] = find(I3);
    II2 = II1-1;
    II2(II2 == 0) = 3;
    
    K1 = sub2ind(size(T),J,II1);
    K2 = sub2ind(size(T),J,II2);
    C1 = vertcat(C1,[T(K1),T(K2)]);

    % Sorting
    C1 = sort(C1,2);
    C1 = unique(C1,'rows');

%     TRI = delaunayTriangulation([Bx1, By1; Cx, Cy],C);
%     figure; triplot(TRI); hold on; line(P(C1',1),P(C1',2));

    %----------------------------------------------------------------------
    % Visualize
    %----------------------------------------------------------------------
%     figure; hold on;
%     triplot(TRI,'color',[1 1 1]*.7);
%     line(x1(C1'),y1(C1'),'color','k');
%     plot(P1(:,1),P1(:,2),'r','Linewidth',2);
%     plot(P2(:,1),P2(:,2),'b','linewidth',2);
%     axis equal tight;
% 
%     figure; hold on;
% %     triplot(T,P(:,1),P(:,2),'b');
%     line(x1(C1'),y1(C1'),'color','b');
%     plot(P1(:,1),P1(:,2),'k','Linewidth',2);
%     plot(P2(:,1),P2(:,2),'r','linewidth',2);
%     axis equal tight;


%     I = ismember(P1,P(C(:,1),:),'rows');
% 
%     temp = P(C(:,1),:);
%     P11 = P1;
%     P11(I,:) = [];

    %----------------------------------------------------------------------
    % Fill missing connections
    %----------------------------------------------------------------------
%     tic;
    C2 = [];
    for i = 1 : size(P1,1)
        if ~any(ismember(C1(:,1),i))
            II1 = find(C1(:,1) < i,1,'last');
            II2 = find(C1(:,1) > i,1,'first');

            if isempty(II1) || isempty(II2)
                continue;
            end

            n1 = C1(II1,1);
            n2 = C1(II2,1);

            II1 = find(C1(:,1) == n1);
            II2 = find(C1(:,1) == n2);

            if isempty(II1) || isempty(II2)
                error('Error in filling missing R2C lines.');
            end

            j = intersect(C1(II1,2),C1(II2,2));
            if isempty(j)
                error('Error in finding common node of R2C lines.');
            end

            % Check intersection with boundary and channel
            I = n1:n2;
            [xi,yi] = polyxpoly(x1([i,j]),y1([i,j]),P1(:,1),P1(:,2));
            if length(xi) > 1
                continue;
            end

%             [xi,yi] = polyxpoly(x1([i,j]),y1([i,j]),P2(:,1),P2(:,2));
%             if length(xi) > 1
%                 continue;
%             end

            C2(end+1,:) = [i, j];

        end
    end

    C3 = [];
    for i = 1 : size(P1,1)
        if nnz(ismember(C1(:,1),i)) > 1
            II1 = find(C1(:,1) == i);
            n1 = min(C1(II1,2));
            n2 = max(C1(II1,2));
            
            if n2 - n1 < 2
                continue;
            end

            J = n1 : n2;
            for j = n1:n2
                % Check intersection with boundary and channel
%                 [xi,yi] = polyxpoly(x1([i,j]),y1([i,j]),P1(:,1),P1(:,2));
%                 if length(xi) > 1
%                     continue;
%                 end

                [xi,yi] = polyxpoly(x1([i,j]),y1([i,j]),P2(:,1),P2(:,2));
                if length(xi) > 1
                    continue;
                end

                C3 = vertcat(C3,[i,j]);

            end
            
            

        end
    end
%     toc;

    %----------------------------------------------------------------------
    % Set global constraint
    %----------------------------------------------------------------------
    C = [C1;C2;C3];
    C = sort(C,2);
    C(:,2) = -1*C(:,2);
    C = unique(C,'rows');
    C(:,2) = -1*C(:,2);

    %----------------------------------------------------------------------
    % Visualize
    %----------------------------------------------------------------------
%     figure; hold on;
%     line(x1(C2'),y1(C2'),'color',[1 1 1]*0.7);
%     line(x1(C3'),y1(C3'),'color',[1 1 1]*0.7);
%     line(x1(C1'),y1(C1'),'color',[1 1 1]*0.0);
%     plot(P1(:,1),P1(:,2),'r','Linewidth',2);
%     plot(P2(:,1),P2(:,2),'b','linewidth',2);
%     axis equal tight;

    %% ----------------------------------------------------------------------
    % Smooth topography along connection line
    %----------------------------------------------------------------------
    PIs = cell(size(C,1),1);
    DataSets = cell(size(C,1),1);
    if nConnectingLines >= size(C,1)
        IC = 1:size(C,1);
    else
        IC = linspace(1,size(C,1),nConnectingLines+1);
        IC = IC(1:end-1);
%         IC = linspace(1,size(C,1),nConnectingLines);
        IC = round(IC);
    end
    

    x1 = P(:,1);
    y1 = P(:,2);

    for i = IC
        % Set connection line
        x2 = linspace(x1(C(i,1)),x1(C(i,2)),1e2);
        y2 = linspace(y1(C(i,1)),y1(C(i,2)),1e2);
        
        % Retreive cross-section topography
        r2 = sqrt( (x2 - x2(1)).^2 + (y2 - y2(1)).^2);
        z2 = -xyzFun(x2',y2');
        z2 = z2';
        
        % Apply cubic spline smoothing until there is no upslope
        RMSE = 0.1;
        while 1
            PI = ComputePathCurvature([r2(:), z2(:)],[r2([1 end])' z2([1 end])'],RMSE,'warning','off');
            PI.RMSE = RMSE;
            if max(PI.dy(PI.p)) < 0
                break;
            end
            if z2(1) <= z2(end)
                break;
            end

            RMSE = RMSE+0.1;
        end
        PIs{i} = PI;
        DataSets{i}.x = x2;
        DataSets{i}.y = y2;
    end

    % Remove empty sets
    I = cellfun(@(x) isempty(x),PIs);
    PIs(I) = [];
    DataSets(I) = [];

    % Set points to be used in interpolation
    x1 = [];
    y1 = [];
    z1 = [];
    for i = 1 : length(PIs)
        PI = PIs{i};
        x1 = [x1; DataSets{i}.x(:)];
        y1 = [y1; DataSets{i}.y(:)];
        z1 = [z1; PI.sy(PI.p)];
    end

    x1 = [x1; Cx; Bx];
    y1 = [y1; Cy; By];
    z1 = [z1; Cz; Bz];

    [~,I] = unique([x1,y1],'rows');
    x1 = x1(I);
    y1 = y1(I);
    z1 = z1(I);

    % Interpolate
    fInterp = scatteredInterpolant(x1,y1,z1,'linear','none');
    Znew1 = fInterp(X,Y);

    % Assign new elevation data within basin mask
    Znew(BasinMask) = Znew1(BasinMask);


    % Visualize
%     figure;
%     for i = 1 : length(PIs)
%         subplot(1,2,1); hold off;
%         plot(P1(:,1),P1(:,2),'k','Linewidth',2);
%         hold on;
%         plot(P2(:,1),P2(:,2),'r','linewidth',2);
%         plot(DataSets{i}.x,DataSets{i}.y,'b');
%         axis equal tight;
% 
%         subplot(1,2,2); hold off
%         PI = PIs{i};
%         plot(PI.x1(PI.p),PI.y1(PI.p));
%         hold on;
%         plot(PI.sx(PI.p),PI.sy(PI.p));
% 
%         drawnow;
%     end
% 
%     figure; hold on;
%     plot(P1(:,1),P1(:,2),'k','Linewidth',2);
%     plot(P2(:,1),P2(:,2),'r','linewidth',2);
%     scatter3(x1,y1,z1,1,z1,'filled');
%     axis equal tight;


    
end

Znew = -Znew;
Znew(isnan(Znew)) = Z(isnan(Znew));

app.xyzFun_new = app.xyzFun;
app.xyzFun_new.Values = Znew;
app.DEMDropDown.Enable = 'on';

PlotR2D2Result(app);
