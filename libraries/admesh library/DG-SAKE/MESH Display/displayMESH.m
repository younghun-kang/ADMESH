function displayMESH(MESH,RO_Stations,arrows)
%==========================================================================
% Display multi-dimensional MESH framework and additional attributes
%==========================================================================
fprintf('Plotting mesh...\n\n');
%--------------------------------------------------------------------------
% Open Figure
%--------------------------------------------------------------------------
fig = figure('color','w'); hold on

MESH.Points(:,3) = -MESH.Points(:,3);
%--------------------------------------------------------------------------
% Compute colormap index for mesh
%--------------------------------------------------------------------------
zmin = min(MESH.Points(:,3)); % Compute min elevation
zmax = max(MESH.Points(:,3)); % Compute max elevation

% Normalize and multiply by 256 (number of rows in colormap)
C = round(256*(MESH.Points(:,3) - zmin) / ( zmax - zmin ));
C(C>256) = 256;  C(C <= 0) = 1; % Cap

% Get landcover colormap
cmap = flipud(landcolor(256));

%--------------------------------------------------------------------------
% Plot MESH
%--------------------------------------------------------------------------
patch(...
    'faces',MESH.ConnectivityList,...
    'vertices',MESH.Points(:,[1 2]),...
    'facevertexcdata',cmap(C,:),...
    'facecolor','interp',...
    'edgecolor','k',...
    'edgealpha',.2);

%--------------------------------------------------------------------------
% Plot derived cascade network and channels
%--------------------------------------------------------------------------
if isfield(MESH,'channelNetwork')
    
    % Get default line colors
    cax = get(gca,'colororder');

    % Plot cascades
    for k = 1:length(MESH.channelNetwork)
        
        % Get coordinates
        x = MESH.Points(MESH.channelNetwork{k}(:,1),1);
        y = MESH.Points(MESH.channelNetwork{k}(:,1),2);
        %z = MESH.Points(MESH.channelNetwork{k}(:,1),3);
        
        % Compute elevation offset
        %offset = z*.001;
        
        % Plot line
        %line(x,y,z + offset,'color',cax(1,:),'linestyle','-','linewidth',2);
        line(x,y,'color','r','linestyle','-','linewidth',2);
        
        % Plot beginning and end nodes
        %plot3(...
        %    x([1,end]),y([1,end]),z([1,end]) + offset([1,end]),...
        %    'color',cax(1,:),'linestyle','none','marker','o',...
        %    'markerfacecolor',cax(6,:),'markersize',4)
        
    end
    
end

%--------------------------------------------------------------------------
% Plot channel network in mesh
%--------------------------------------------------------------------------
for k = 1:length(MESH.Constraints)
    if MESH.Constraints(k).num == 18 || MESH.Constraints(k).num == 6
        % Get coordinates
        x = MESH.Points(MESH.Constraints(k).nodeStr,1);
        y = MESH.Points(MESH.Constraints(k).nodeStr,2);
        %z = MESH.Points(MESH.Constraints(k).nodes,3);
        
        % Compute elevation offset
        %offset = z*.001;
        
        % Plot line
        %line(x,y,z + offset,'color',cax(3,:),'linestyle','--','linewidth',2);
        if all(MESH.Constraints(k).data(:,3) ~= 0)
        line(x,y,'color','b','linestyle','-','linewidth',2);
        else
        line(x,y,'color','r','linestyle','-','linewidth',2);
        end
        % Plot beginning and end nodes
%         plot3(...
%             x([1,end]),y([1,end]),z([1,end]) + offset([1,end]),...
%             'color',cax(1,:),'linestyle','none','marker','o',...
%             'markerfacecolor',cax(6,:),'markersize',4)
    end
end


% Compute unit normal vector to each triangle
if arrows == 1
    [ts,~,~] = triSlopeVector(MESH);
    
    % Compute centroid
    ic = (...
        MESH.Points(MESH.ConnectivityList(:,1),:) + ...
        MESH.Points(MESH.ConnectivityList(:,2),:) + ...
        MESH.Points(MESH.ConnectivityList(:,3),:))/3;
    
    % Compute elevation offset
    offset = ic(:,3)*.001;

    for i = 1 : size(MESH.ConnectivityList,1)
    [~,A(i)] = boundary(MESH.Points(MESH.ConnectivityList(i,:),1),MESH.Points(MESH.ConnectivityList(i,:),2));
    end

    ts1 = sqrt(sum(ts.^2,2));
    ts = ts./ts1.*sqrt(A(:));
    
    % Display flow direction
    %quiver3(ic(:,1),ic(:,2),ic(:,3)+offset,ts(:,1),ts(:,2),ts(:,3),.5,'color','k')
    quiver(ic(:,1),ic(:,2),ts(:,1),ts(:,2),1,'color','k')
    
end

%==========================================================================
% Plot run off stations
%==========================================================================
if ~isempty(RO_Stations)
for k = 1:size(RO_Stations.id,1)
    
    % Assign coordinates
    x = RO_Stations.xs(k);
    y = RO_Stations.ys(k);
    %z = xyzFun(x,y);
    
    % Compute elevation offset
    %offset = z*.1;
    
    % Plot station marker
    plot(x,y,'ko','markerfacecolor','w','markersize',10);
    
    % Plot station number
    text(x,y,num2str(RO_Stations.id(k)),'color','k','fontweight','bold','fontsize',18);
    
end
end



% for k = 1:size(MESH.ConnectivityList,1)
%     
%     % Plot main element
%     pp = patch(...
%         'faces',MESH.ConnectivityList(k,:),...
%         'vertices',MESH.Points,...
%         'facecolor','r',...
%         'edgecolor','k',...
%         'edgealpha',.1);
%     
%     if ~isempty(MESH.flowMap{k})
%         % Plot elements contributing flow
%         pq = patch(...
%             'faces',MESH.ConnectivityList(MESH.flowMap{k},:),...
%             'vertices',MESH.Points,...
%             'facecolor','g',...
%             'edgecolor','k',...
%             'edgealpha',.6);
%     end
%     pause
%     
%     delete([pp pq])
%     
%     
% end

%==========================================================================
% Plot rain guage stations
%==========================================================================
% L = load('RainGuageStationLocation.pre');
% 
% for k = 1:size(L,1)
%     
%     % Assign coordinates
%     x = L(k,2);
%     y = L(k,3);
%     z = xyzFun(x,y);
%     
%     % Compute elevation offset
%     offset = z*.1;
%     
%     % Plot station marker
%     plot3(x,y,z + offset,'ko','markerfacecolor',cax(2,:),'markersize',10);
%     
%     % Plot station number
%     text(x,y,z + offset,num2str(L(k,1)),'color','k','fontweight','bold','fontsize',18);
%     
% end

axis equal;
return;
%--------------------------------------------------------------------------
% Set axis properties
%--------------------------------------------------------------------------
colormap(cmap)
daspect([1 1 1/13])
xmin = min(MESH.Points(:,1)); xmax = max(MESH.Points(:,1));
ymin = min(MESH.Points(:,2)); ymax = max(MESH.Points(:,2));
zmin = min(MESH.Points(:,3)); zmax = max(MESH.Points(:,3));
offsetx = abs(xmin - xmax)*.1;
offsety = abs(ymin - ymax)*.1;
axis([(xmin-offsetx) (xmax+offsetx) (ymin-offsety) (ymax+offsety)])

cbh = colorbar;
cbh.YTick = linspace(0,.99,10);
set(cbh,'yticklabel',arrayfun(@(x) num2str(x,'%.0f'),round(linspace(zmin, zmax,10)),'uni',false));
set(get(cbh,'Title'),'String','m')
xlabel('X (m)'), ylabel('Y (m)'), zlabel('Z (m)')
box on
set(gca,'fontsize',36,'fontname','Times')
plot_border
view(3)
maxfig(fig,1); drawnow

end
