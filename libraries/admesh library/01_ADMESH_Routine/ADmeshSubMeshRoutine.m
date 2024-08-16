function ADmeshSubMeshRoutine(varargin)
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
gui = % guidata(varargin{3});

%--------------------------------------------------------------------------
% Check to make sure user an edge structure has been uploaded
%--------------------------------------------------------------------------
quit = CheckSubMeshUserInput(gui.Window);

% Stop ADMESH routine if something's not right
if quit == 1 
    return
end

% Turn off colormap
SetContourStatus(gui.Window,nan,'off'); drawnow

%--------------------------------------------------------------------------
% Assign ADmesh parameters
%--------------------------------------------------------------------------
[res,hmax,hmin,K,R,s,C,T,g] = DealVariables(gui.Window);

%--------------------------------------------------------------------------
% Save ADMESH settings for appending to file
%--------------------------------------------------------------------------
Settings = SaveSettings(gui.Window); 

%--------------------------------------------------------------------------
% Begin ADmesh
%--------------------------------------------------------------------------
t=tic; % Start Timer for ADmesh

%--------------------------------------------------------------------------
% Create structured background mesh
%--------------------------------------------------------------------------
[X,Y,delta] = CreateBackgroundGrid(gui.subPTS,hmax,hmin,res,gui.sb);

%--------------------------------------------------------------------------
% Calculate distance function
%--------------------------------------------------------------------------
[D,gradD] = SignedDistanceFunction(gui.subPTS,X,Y,delta,hmax,gui.sb);

% figure
% contourf(X,Y,D)
% colorbar
% return

%--------------------------------------------------------------------------
% Compute boundary curvature
%--------------------------------------------------------------------------
h_curve = CurvatureFunction(D,gradD,X,Y,K,g,hmax,hmin,Settings,gui.sb);

d = (bwdist(D < 0) - bwdist(D > 0))*hmin;

ix = find(abs(d) <= 4*hmin);

h_curve(ix(h_curve(ix) < 900)) = 900;
%--------------------------------------------------------------------------
% Compute local feature size
%--------------------------------------------------------------------------
h_lfs = MedialAxisFunction(X,Y,D,gradD,R,hmin,hmax,Settings,gui.sb);

%--------------------------------------------------------------------------
% Compare initial conditions and save
%--------------------------------------------------------------------------
h_ic = min(h_curve, h_lfs); clear h_curve h_lfs

%--------------------------------------------------------------------------
% Interpolate bathymetry to background grid if needed.
%--------------------------------------------------------------------------
Z = CreateElevationGrid(X,Y,gui.xyzFun,Settings,gui.sb);

%--------------------------------------------------------------------------
% Compute bathymetry
%--------------------------------------------------------------------------
h_bathy = BathymetryFunction(X,Y,Z,s,hmin,hmax,delta,Settings,gui.sb);
h_bathy(h_bathy < 700) = 700;
%--------------------------------------------------------------------------
% Compare initial conditions and save
%--------------------------------------------------------------------------
h_ic = min(h_ic, h_bathy); clear h_bathy

%--------------------------------------------------------------------------
% Compute dominate tide function
%--------------------------------------------------------------------------
h_tide = Dominate_tide(T,C,Z,size(X),hmax,hmin,Settings); clear Z

%--------------------------------------------------------------------------
% Compare initial conditions and save
%--------------------------------------------------------------------------
h_ic = min(h_ic, h_tide); clear h_tide

%--------------------------------------------------------------------------
% Enforce boundary conditions
%--------------------------------------------------------------------------
h_ic = EnforceBoundaryConditions(h_ic,X,Y,D,gui.subPTS,hmax,hmin,gui.sb);

%--------------------------------------------------------------------------
% Compute mesh size function
%--------------------------------------------------------------------------
% gui.sb.setText('Computing Mesh Size Function...')
h = MeshSizeIterativeSolver(h_ic,D,hmax,hmin,g,delta);
% gui.sb.setText('Preparing to start mesh generation...')

%--------------------------------------------------------------------------
% Create an interpolant surface for the distance and mesh size function 
%--------------------------------------------------------------------------
phi.f       = griddedInterpolant(X',Y',D'); 
phi.fx      = griddedInterpolant(X',Y',gradD.x'); 
phi.fy      = griddedInterpolant(X',Y',gradD.y'); clear gradD D  
MeshFun     = griddedInterpolant(X',Y',h'); clear X Y 
hmin        = min(h(:)); clear h

%--------------------------------------------------------------------------
% Generate mesh
%--------------------------------------------------------------------------
submesh = distsubmesh2d(gui.subPTS,phi,MeshFun,gui.xyzFun,hmin,Settings,gui.sb,gui.ViewAxes);

clear DistFun DistGxFun DistGyFun MeshFun

%--------------------------------------------------------------------------
% Plot Final Results
%--------------------------------------------------------------------------
% gui.sb.setText('Displaying final mesh...')
PlotSubMesh(gui,submesh)

%--------------------------------------------------------------------------
% Convert CPU time to hours minutes seconds time string
%--------------------------------------------------------------------------
time_string = seconds2HrMinSec(toc(t));

%--------------------------------------------------------------------------
% Update status bar
%--------------------------------------------------------------------------
% gui.sb.setText(['Run time: ' , time_string ])

% Turn on coordinate display. 
set(gui.Window,'WindowButtonMotionFcn' , @CoordDisplay)

%--------------------------------------------------------------------------
% Ask user how it looks?
%--------------------------------------------------------------------------
msg = 'Would you like to keep the sub mesh or run ADMESH again?';
choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'Continue','Re-Try'},'DefaultOption',1,'Icon','Warning');

switch choice
    
    case 'Continue'
        
        % Append sub mesh to full mesh
        
    case 'Re-Try'
        
        % Re-plot sub boundary
        
        % Delete current plot
        h = findobj(gui.ViewAxes,'tag', 'sub mesh'); if ~isempty(h); delete(h); end
        % Plot polygon
        for k = 1:length(gui.subPTS.Poly)
            
            hold on
            bH = plot(gui.subPTS.Poly(k).x,gui.subPTS.Poly(k).y,'linewidth',2,'Color',[.7 .5 0]);
            set(bH ,'Tag','sub domain');
            
        end
        
        % Plot subdomain constraint
        for k = 1:length(gui.subPTS.Constraints)
            
            hold on
            p = gui.subPTS.Constraints(k).xy;
            bH = plot(p(:,1),p(:,2),'linewidth',2,'Color','r');
            set(bH ,'Tag','sub domain');
            
        end
        
end

return

%--------------------------------------------------------------------------
% Update guidata
%--------------------------------------------------------------------------
% guiH = % guidata(guiFig);
% 
% guiH.mesh = mesh; 
% 
% guidata(guiFig,guiH);

%--------------------------------------------------------------------------
% Write submesh to kml
%--------------------------------------------------------------------------
MESH = submesh;

MESH.cpplon = gui.MESH.cpplon;
MESH.cpplat = gui.MESH.cpplat;

%--------------------------------------------------------------------------
% Check for variables
%--------------------------------------------------------------------------
if isempty(MESH) % User has not run ADmesh yet
    msg = 'No mesh to save....';
    uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');
    return
end

%--------------------------------------------------------------------------
% Prepare mesh for kml write out. 
%--------------------------------------------------------------------------

% Convert coordinates
if isfield(MESH,'cpplon')
    MESH = Meters2Geo(MESH);
else
    msg = ['In order to write a kml file the coordinates ',...
        'of the domain needed to be geographic originally.'];
    uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');
    return
end

% Ask user for file name
% gui.sb.setText('Save file as...')
[file,path] = uiputfile('*.kml','Save As');

% If user cancels
if isequal(file,0) || isequal(path,0)
    % gui.sb.setText('Ready')
    return
end

%--------------------------------------------------------------------------
% Break down each element into line segments
%--------------------------------------------------------------------------
% Assign nodes and connectivity list
p = MESH.Points;
t = MESH.ConnectivityList;

% Get triangulation representation
trep = triangulation(t,p);

% Get edges, sort left to right
E = sort(edges(trep),2);

% Initialize cell for holding boundary conditions
BC = cell(numel(MESH.Constraints),1);

% Remove boundary conditions from edges matrix
for k = 1:numel(MESH.Constraints)
    
    % Need to handle external barriers differently
    if any(MESH.Constraints(k).num == [4 5 24 25])
        
        % 1st Node string
        nStr1 = MESH.Constraints(k).nodeStr(:,1);
        % 2nd Node string
        nStr2 = MESH.Constraints(k).nodeStr(:,2);
        
        % Compete edge list
        ns = [[nStr1(1:end-1),nStr1(2:end)];[nStr2(1:end-1),nStr2(2:end)]];
        
        % Append edge type
        type = ones(size(ns,1),1)*MESH.Constraints(k).num;
        
        % Constraint (append type)
        BC{k} = [sort(ns,2) type];

    else
        
        % Node string
        ns = MESH.Constraints(k).nodeStr;
        
        % Append edge type
        type = ones(size(ns,1)-1,1)*MESH.Constraints(k).num;
        
        % Constraint (append type)
        BC{k} = [sort([ns(1:end-1), ns(2:end)],2) type];
        
    end
    
    E = setdiff(E,BC{k}(:,[1 2]),'rows'); % Remove common points
    
end

% For regular element edges, append nan
E = [E nan(size(E,1),1)];

% Concatenate elements edges
E = [cell2mat(BC); E];

%--------------------------------------------------------------------------
% Determine if file should be broken up
%--------------------------------------------------------------------------
% var = whos('E');
% 
% if var.bytes > 1000 % Break up
%     
%     
%     
%     
%     
% end


%--------------------------------------------------------------------------
% Prepare to write
%--------------------------------------------------------------------------

% Initialize kml header
Header = {...
    '<?xml version="1.0" encoding="UTF-8"?>',...
    '<kml xmlns="http://www.opengis.net/kml/2.2">',...
    '   <Document>',...
    '   <Style id="yellowStyle">',...
    '       <LineStyle>',...
    '           <color>ff00ffff</color>',...
    '           <width>1</width>',...
    '       </LineStyle>',...
    '   </Style>',...
    '   <Style id="redStyle">',...
    '       <LineStyle>',...
    '           <color>ff0000ff</color>',...
    '           <width>2</width>',...
    '       </LineStyle>',...
    '   </Style>',...
    '   <Style id="blueStyle">',...
    '       <LineStyle>',...
    '           <color>ffff0000</color>',...
    '           <width>2</width>',...
    '       </LineStyle>',...
    '   </Style>',... 
    '   <Style id="greenStyle">',...
    '       <LineStyle>',...
    '           <color>ff00ff00</color>',...
    '           <width>2</width>',...
    '       </LineStyle>',...
    '   </Style>',... 
    '   <Style id="blackStyle">',...
    '       <LineStyle>',...
    '           <color>ff000000</color>',...
    '           <width>2</width>',...
    '       </LineStyle>',...
    '   </Style>',...
    '   <Style id="brownStyle">',...
    '       <LineStyle>',...
    '           <color>ff336699</color>',...
    '           <width>2</width>',...
    '       </LineStyle>',...
    '   </Style>'}';

% Initialize beginning & ending block used for each path
bBlock = {...
    '   <Placemark>',...
    '   <name>Element</name>',...
    '   <styleUrl>#yellowStyle</styleUrl>',...
    '   <LineString>',...
    '       <tessellate>1</tessellate>',...
    '       <coordinates>'}';

eBlock = {...
    '       </coordinates>',...
    '   </LineString>',...
    '   </Placemark>'}';

% Initialize kml footer
Footer = {...
    '   </Document>',...
    '</kml>'}';

% Initialze cell for fprintf 
nedges     = size(E,1);
nHeaders   = size(Header,1);
nLines     = 11*nedges;
nFooters   = 2;

textOut = cell(nHeaders + nLines + nFooters,1);

% gui.sb.setText('Preparing mesh for kml write out...')

% Store header
textOut(1:nHeaders) = Header;

% Initialize start and end index for edge write out
s = nHeaders+1;
e = s + 10;

% Determine appropraite frequency for updating status bar
freq = round(nedges/100);

gui.sb.ProgressBar.setVisible(true)
set(gui.sb.ProgressBar, 'Minimum',1, 'Maximum',nedges, 'Value',1)

% Loop over each edge and prep for write out
for k = 1:nedges
    
    % Assign (x,y,z) coordinated for edge
    x = [p(E(k,1),1); p(E(k,2),1)];
    y = [p(E(k,1),2); p(E(k,2),2)];
    z = [p(E(k,1),3); p(E(k,2),3)];
    
    % What color are we using?
    if isnan(E(k,3))
        bBlock(3) = {'   <styleUrl>#yellowStyle</styleUrl>'};
    elseif E(k,3) == -1
        bBlock(3) = {'   <styleUrl>#blueStyle</styleUrl>'};
    elseif any(E(k,3) == [0 2 10 12 20 22 30])
        bBlock(3) = {'   <styleUrl>#blackStyle</styleUrl>'}; 
    elseif any(E(k,3) == [1 11 21])
        bBlock(3) = {'   <styleUrl>#greenStyle</styleUrl>'}; 
    elseif any(E(k,3) == [3 13 23 4 5 24 25])
        bBlock(3) = {'   <styleUrl>#brownStyle</styleUrl>'}; 
    elseif E(k,3) == 18
        bBlock(3) = {'   <styleUrl>#redStyle</styleUrl>'};         
    end
    
    textOut(s:e) = [...
        bBlock; ...
        [{sprintf('               %5.10f,%5.10f,%f',x(1),y(1),z(1) )},...
        {sprintf('               %5.10f,%5.10f,%f',x(2),y(2),z(2) )}]';...
        eBlock];
    
    s = e + 1;
    e = s + 10;
    
    if(mod(k,freq ) == 0)
        set(gui.sb.ProgressBar,'Value',k)
    end
    
end

% Append footer
textOut(end-1:end) = Footer;

gui.sb.ProgressBar.setVisible(false)

%--------------------------------------------------------------------------
% Open file and write
%--------------------------------------------------------------------------

% gui.sb.setText('Writing kml file...')

% Open file
fid = fopen([path file],'w');

fprintf(fid,'%s\n',textOut{:});

% Close file
fclose(fid);

% gui.sb.setText('Ready')


end