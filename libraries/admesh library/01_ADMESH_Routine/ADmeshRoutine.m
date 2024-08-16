function ADmeshRoutine(varargin)
% ADvanced Mesh Generator for hydrodynmaic models 1.1.6 (ADMESH)
%
% Variables:
%   Edge Structure Variable:
%    PTS : Data structure containing polygon (x & y) points of boundaries
%          and boundary conditions
%
%   User specified Variables:
%    res    : a factor multiple, either 1,3,or 4, of hmin controlling the
%             grid spacing.
%    hmin   : Minimum element size
%    hmax   : Maximum element size
%    K      : parameter specifying the number of elements per radian of
%             curvature
%    R      : parameter specifying the number of elements that span half
%             a channel/tributary
%    xyz    : gridded interpolant of elevation
%    s      : parameter specifying the density of elements in areas with
%             large changes in elevation
%    C      : parameter specifying the number of elements per tidal
%             wavelength
%    T      : frequency of wave propogation
%    g      : grading limiting
%
%   Function Output Variables
%    X           : x values of background grid
%    Y           : y values of background grid
%    delta       : cellsize (hmin/res)
%    D           : distance function ('-' values = region to be meshed)
%    gradD       : gradient of distance function
%    Z           : bathymetry on background mesh
%    h_curve     : mesh size values determined via curvature
%    h_lfs       : mesh size values determined via local feature size
%    h_bathy     : mesh size values determined via horizontal length scale
%    h_tide      : mesh size values determined via linear wave speed
%    DistFun     : gridded interpolant of distance function
%    DistGxFun   : gridded interpolant of gradient (x-direction) function
%    DistGyFun   : gridded interpolant of gradient (y-direction) function
%    MeshFun     : gridded interpolant of mesh size function
%    h           : mesh resolution determined via solution of PDE
%    mesh        : Data structure containing x and y coordinates of
%                  unstructured nodes, connectivity of mesh and node
%                  strings of boundary conditions.
%    minEQ       : minimum element quality
%    meanEQ      : mean element quality
%--------------------------------- BEGIN CODE -----------------------------

%--------------------------------------------------------------------------
% Get GUI data
%--------------------------------------------------------------------------
app = varargin{1};

% Turn off colormap
SetContourStatus(app,'Off'); drawnow

Compute1DProjection(app);
SetDummyConstraint(app);

%--------------------------------------------------------------------------
% Check to make sure user an edge structure has been uploaded
%--------------------------------------------------------------------------
quit = CheckUserInput(app);

% Stop ADMESH routine if something's not right
if quit == 1
    return
end

%--------------------------------------------------------------------------
% Assign ADmesh parameters
%--------------------------------------------------------------------------
[res,hmax,hmin,K,R,s,C,T,g] = DealVariables(app);

%--------------------------------------------------------------------------
% Save ADMESH settings for appending to file
%--------------------------------------------------------------------------
[Settings,AdvSettings] = SaveSettings(app);

%--------------------------------------------------------------------------
% Begin ADmesh
%--------------------------------------------------------------------------
t=tic; % Start Timer for ADmesh

%--------------------------------------------------------------------------
% Create structured background mesh
%--------------------------------------------------------------------------
[X,Y,delta] = CreateBackgroundGrid(app.PTS,hmax,hmin,res,app.UIFigure);

%--------------------------------------------------------------------------
% Calculate distance function
%--------------------------------------------------------------------------
[D,gradD]   = SignedDistanceFunction(app.PTS,X,Y,delta,hmax,app.UIFigure);

while 1 % Repeat if 1D mesh generation turns on
%--------------------------------------------------------------------------
% Initialize initial condition for mesh size function
%--------------------------------------------------------------------------
h0 = hmax*ones(size(D));

%--------------------------------------------------------------------------
% Compute boundary curvature
%--------------------------------------------------------------------------
h0 = CurvatureFunction(h0,D,gradD,X,Y,K,g,hmax,hmin,Settings,app.UIFigure,app.PTS);

%--------------------------------------------------------------------------
% Compute local feature size
%--------------------------------------------------------------------------
h0 = MedialAxisFunction(h0,X,Y,D,gradD,R,hmin,hmax,Settings,app.UIFigure);

%--------------------------------------------------------------------------
% Interpolate bathymetry to background grid if needed.
%--------------------------------------------------------------------------
Z  = CreateElevationGrid(X,Y,app.xyzFun,Settings,app.UIFigure);

%--------------------------------------------------------------------------
% Compute bathymetry
%--------------------------------------------------------------------------
h0 = BathymetryFunction(h0,X,Y,Z,s,hmin,hmax,delta,Settings,app.UIFigure);

%--------------------------------------------------------------------------
% Compute dominate tide function
%--------------------------------------------------------------------------
h0 = Dominate_tide(h0,T,C,Z,size(X),hmax,hmin,Settings); clear Z

%--------------------------------------------------------------------------
% Enforce boundary conditions
%--------------------------------------------------------------------------
% h0 = EnforceBoundaryConditions(h0,X,Y,D,gui.PTS,hmax,hmin,gui.sb);

%--------------------------------------------------------------------------
% Compute mesh size function
%--------------------------------------------------------------------------
h  = MeshSizeFunction(h0,D,hmax,hmin,g,delta,app.UIFigure); clear h0

% h(D > 0) = nan;
% 
% figure
% contourf(X,Y,h)
% colorbar
% 
% return

%--------------------------------------------------------------------------
% 2021-04-07 YK: save updated target mesh size to generate "second" 1D mesh
%--------------------------------------------------------------------------
if app.DummyConstraint == 0 % Break while loop if not dummy
    
    break;
    
elseif app.DummyConstraint == 1 % Repeat loop with updated constraints if dummay
    
    GenerateMeshOnConstraints(X,Y,h,Settings,app);

end

end

%--------------------------------------------------------------------------
% Create an interpolant surface for the distance and mesh size function
%--------------------------------------------------------------------------
% gui.sb.setText('Preparing to start mesh generation...')
phi.f       = griddedInterpolant(X',Y',D');
phi.fx      = griddedInterpolant(X',Y',gradD.x');
phi.fy      = griddedInterpolant(X',Y',gradD.y');   clear gradD D
MeshFun     = griddedInterpolant(X',Y',h');         clear X Y
hmin        = min(h(:));                            clear h

%--------------------------------------------------------------------------
% Generate mesh
%--------------------------------------------------------------------------
%MESH = distquadmesh2d(PTS,DistFun,DistGxFun,DistGyFun,MeshFun,xyzFun,hmin,Settings);
app.MESH = distmesh2d(app.PTS,phi,MeshFun,app.xyzFun,hmin,Settings,app.UIFigure,app.UIAxes,AdvSettings);
clear phi MeshFun

%--------------------------------------------------------------------------
% Update guidata
%--------------------------------------------------------------------------
% guidata(fig,gui);

%--------------------------------------------------------------------------
% Plot Final Results
%--------------------------------------------------------------------------
msg = 'Displaying final mesh...';
uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');
PlotMesh(app,.1)

%--------------------------------------------------------------------------
% Convert CPU time to hours minutes seconds time string
%--------------------------------------------------------------------------
time_string = seconds2HrMinSec(toc(t));

%--------------------------------------------------------------------------
% Update status bar
%--------------------------------------------------------------------------
app.ResultsBox.Value{end+1} = '';
app.ResultsBox.Value{end+1} = ['Run time: ' , time_string ];

%--------------------------------------------------------------------------
% Update .mat file
%--------------------------------------------------------------------------
if ~isempty(app.FilePath)
    save(app.FilePath,'Settings', '-append' );
end

% Turn on coordinate display.
% set(app.Window,'WindowButtonMotionFcn' , @CoordDisplay)




%--------------------------------------------------------------------------
% Screwing around
%--------------------------------------------------------------------------

ext_plot = 0;

if ext_plot == 1
    
    figure('units','normalized','pos',[.23,.15,.45,.7],'color','w'); hold on    % Initialize figure
    
    patchinfo.Vertices              = app.MESH.Points;
    patchinfo.Faces                 = app.MESH.ConnectivityList;
    patchinfo.EdgeColor             = 'b';
    patchinfo.FaceColor             = 'none';
    patchinfo.linewidth             = 1.0;
    patch(patchinfo);
    xlabel('X (m)','FontWeight','bold','FontSize',10)                             % X-Label
    ylabel('Y (m)','FontWeight','bold','FontSize',10)                             % Y-Label
    daspect([1 1 1])                            % Data aspect ratio
    
    xmin = min(app.MESH.Points(:,1));
    xmax = max(app.MESH.Points(:,1));
    ymin = min(app.MESH.Points(:,2));
    ymax = max(app.MESH.Points(:,2));
    
    XMIN = xmin - (xmax-xmin)*.05;
    XMAX = xmax + (xmax-xmin)*.05;
    YMIN = ymin - (ymax-ymin)*.05;
    YMAX = ymax + (ymax-ymin)*.05;
    
    axis([XMIN XMAX YMIN YMAX])
    
    hold on
    
    % Plot constraints
    for k = 1:numel(app.MESH.Constraints)
        
        %    Determine line color
        if app.MESH.Constraints(k).num == -1 % Open
            
            %         nStr = gui.MESH.Constraints(k).nodeStr;
            %
            %         patchinfo.Faces             = [ nStr(1:end-1),nStr(2:end) ];
            %         patchinfo.EdgeColor         = 'b';
            %         patchinfo.linewidth         = 1.5;
            %         h                           = patch(patchinfo);
            %         set(h,'tag','Mesh')
            
        elseif any(app.MESH.Constraints(k).num == [0 2 10 12 20 22 30]); % External Boundary
            
            %         nStr = gui.MESH.Constraints(k).nodeStr;
            %
            %         patchinfo.Faces             = [ nStr(1:end-1),nStr(2:end) ];
            %         patchinfo.EdgeColor         = 'k';
            %         patchinfo.linewidth         = 1.5;
            %         h                           = patch(patchinfo);
            %         set(h,'tag','Mesh')
            
        elseif any(app.MESH.Constraints(k).num == [1 11 21]); % Internal Boundary
            
            %         nStr = gui.MESH.Constraints(k).nodeStr;
            %
            %         patchinfo.Faces             = [ nStr(1:end-1),nStr(2:end) ];
            %         patchinfo.EdgeColor         = 'g';
            %         patchinfo.linewidth         = 1.5;
            %         h                           = patch(patchinfo);
            %         set(h,'tag','Mesh')
            
        elseif any(app.MESH.Constraints(k).num == [3 13 23]); % External Constraints
            %
            %         nStr = gui.MESH.Constraints(k).nodeStr;
            %
            %         patchinfo.Faces             = [ nStr(1:end-1),nStr(2:end) ];
            %         patchinfo.EdgeColor         = [.7 .5 0];
            %         patchinfo.linewidth         = 1.5;
            %         h                           = patch(patchinfo);
            %         set(h,'tag','Mesh')
            
        elseif any(app.MESH.Constraints(k).num == 18); % Channel Constraints
            
            nStr = app.MESH.Constraints(k).nodeStr;
            
            patchinfo.Faces             = [ nStr(1:end-1),nStr(2:end) ];
            patchinfo.EdgeColor         = 'r';
            patchinfo.linewidth         = 1.5;
            h                           = patch(patchinfo);
            set(h,'tag','Mesh')
            
        elseif any(app.MESH.Constraints(k).num == [4 5 24 25]); % Internal Constraints
            
            %        % Plot first side
            %         nStr = gui.MESH.Constraints(k).nodeStr(:,1);
            %
            %         patchinfo.Faces             = [ nStr(1:end-1),nStr(2:end) ];
            %         patchinfo.EdgeColor         = [.7 .5 0];
            %         patchinfo.linewidth         = 1.5;
            %         h                           = patch(patchinfo);
            %         set(h,'tag','Mesh')
            %
            %         %Plot second side
            %         nStr = gui.MESH.Constraints(k).nodeStr(:,2);
            %
            %         patchinfo.Faces             = [ nStr(1:end-1),nStr(2:end) ];
            %         patchinfo.EdgeColor         = [.7 .5 0];
            %         patchinfo.linewidth         = 1.5;
            %         h                           = patch(patchinfo);
            %         set(h,'tag','Mesh')
            %
        end
        
    end
    
    % Modify current axes
    set(gca,'Ticklength',[0 0],'box', 'on','xcolor' ,'k','ycolor' ,'k')
    
    % Get axes limits
    xlim = get(gca,'XLim'); ylim = get(gca,'YLim');
    
    % Create borders
    
    % Compute border widths
    bwy = diff(ylim)/100;  bwx = bwy;
    
    % transparent borders
    patch([xlim(1)+bwx,xlim(2)-bwx,xlim(2)-bwx,xlim(1)+bwx],ylim(1) + bwy*[0,0,1,1],'k','FaceColor','none','clipping','off')
    patch([xlim(1)+bwx,xlim(2)-bwx,xlim(2)-bwx,xlim(1)+bwx],ylim(2) - bwy*[0,0,1,1],'k','FaceColor','none','clipping','off')
    patch(xlim(1) + bwx*[0,0,1,1],[ylim(1)+bwy,ylim(2)-bwy,ylim(2)-bwy,ylim(1)+bwy],'k','FaceColor','none','clipping','off')
    patch(xlim(2) - bwx*[0,0,1,1],[ylim(1)+bwy,ylim(2)-bwy,ylim(2)-bwy,ylim(1)+bwy],'k','FaceColor','none','clipping','off')
    
    % Number of bars
    nB = 5; xspace = linspace(xlim(1),xlim(2),nB*2);
    
    for k = 2:2:nB*2
        
        xp = [xspace(k-1) xspace(k-1) xspace(k) xspace(k)];
        patch(xp,ylim(1) + bwy*[0,1,1,0],'k','FaceColor','k','clipping','off')
        patch(xp,ylim(2) - bwy*[0,1,1,0],'k','FaceColor','k','clipping','off')
        
    end
    
    % Number of bars
    nB = 5;
    yspace = linspace(ylim(1),ylim(2),nB*2);
    
    for k = 2:2:nB*2
        
        yp = [yspace(k-1) yspace(k) yspace(k) yspace(k-1)];
        patch(xlim(1) + bwx*[0,0,1,1],yp,'k','FaceColor','k','clipping','off')
        patch(xlim(2) - bwx*[0,0,1,1],yp,'k','FaceColor','k','clipping','off')
        
    end
    
end

end
