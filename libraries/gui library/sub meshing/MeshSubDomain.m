function MeshSubDomain(varargin)
% ExtractSubDomain - Extracts subdomain from PTS data structure
%
% Syntax:  ExtractSubDomain(varargin)
%
% Inputs:
%    None
%
% Outputs:
%    None
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% Last revision: 29-April-2014

%---------------------------------------------------------------------
% Begin Code
%---------------------------------------------------------------------

%---------------------------------------------------------------------
% Get GUI Data
%---------------------------------------------------------------------
app = varargin{1};

% Assign variables
PTS     = app.PTS;
MESH    = app.MESH;

%---------------------------------------------------------------------
% Check if any data is loaded
%---------------------------------------------------------------------
pH = app.UIAxes; % Plot Window handle

% Catch no data loaded
if isempty(MESH) && isempty(PTS)
    
    errordlg('No file is currently loaded.','Error');
    return
    
end

% Catch no mesh
if isempty(MESH) && ~isempty(PTS)
    
    errordlg('You must first generate a mesh.','Error');
    return
    
end

% Catch no PTS
if ~isempty(MESH) && isempty(PTS)
    
    if ~isempty(MESH.Constraints)
        PTS = Constraints2PTS(MESH);
    else
        PTS = triangulation2PTS(MESH);
    end
    
    app.PTS = PTS;
    %% guidata(gui.Window,gui)
    
end

% Turn off colormap
SetContourStatus(app,'Off'); drawnow

% If user has a mesh loaded but no PTS strucutre, make one
if ~isempty(MESH) && ~isempty(PTS)
    
    msg = [...
        'ADMESH will be using the current edge structure to ',...
        'define the external and internal boundaries. If there ',...
        'is a more refined boundary you would like to use then ',...
        'select from the options below to choose an alternative ',...
        'boundary to use.'];
    choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'Continue','Load Alternative Boundary','Quit'},'DefaultOption',1,'Icon','Warning');
    drawnow; pause(.005);
    
    switch choice
        
        
        case 'Continue'
            
            % Do nothing, continue
            
        case 'Load Alternative Boundary'
            
            msg = ['You have several options to choose from. ',...
                '(1) Load an existing edge structure (.mat) file. ',...
                '(2) Import boundary from another mesh (.14) file. ',...
                '(3) Import a polygon from a google earth kml file.'];
            choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
                'Options',{'(.mat) File','(.14) File','(.kml) File'},'DefaultOption',1,'Icon','Warning');
            switch choice
                
                case '(.mat) File'
                    
                    fileType = '.mat';
                    
                case '(.14) File'
                    
                    fileType = '.14';
                    
                case '(.kml) File'
                    
                    fileType = '.kml';
                    
                case ''
                    
                    return
  
            end
            
            [app.PTS,status,xyzIncluded,elev] = ImportTemporaryFile(fileType,MESH);
            
            if ~status; % gui.sb.setText('Continue?'); return; end
            
            if xyzIncluded; app.xyzFun = elev; clear elev; end

            % update gui data
            % guidata(fig,gui);
            
            PTS     = app.PTS;            
            
        case 'Quit'

            return
            
        case '' % User exit box

            return

    end
        
end

% Set all toolbar tools off
%set(gui.uiToolbarTool([1,2,3]),'State','off')

%---------------------------------------------------------------------
% Plot Actual boundary
%---------------------------------------------------------------------
axes(pH); hold on
bH = line(...
    'Parent',pH,...
    'Xdata',PTS.Poly(1).x,...
    'Ydata',PTS.Poly(1).y,...
    'LineWidth',2,...
    'Color','r');
drawnow

%---------------------------------------------------------------------
% Select mesh domain
%---------------------------------------------------------------------
% gui.sb.setText('Finding elements in selected region...')

% Current mesh 
[xnodes,ynodes,t] = deal(MESH.Points(:,1),MESH.Points(:,2),MESH.ConnectivityList);

% Ask user to draw polygon around region
bpoly = DrawBoundingPolygon(gui.Window,pH,'r'); %assignin('base','bpoly',bpoly)

% Find vertices in polygon
verts = find(PointInPolygon(xnodes,ynodes,bpoly(:,1), bpoly(:,2)));

% Get triangulation representation of current mesh
trep = triangulation(t,xnodes,ynodes);

% Find all elements attached to vertices
ti = vertexAttachments(trep,verts); ti = [ti{:}]';

% Remove elements containing vertices that do not fall within the selected
% polygon.
ti = ti(sum(ismember(t(ti,:),verts),2) == 3);

% gui.sb.setText('Highlighting elements in selected region...')

%---------------------------------------------------------------------
% Highlight nodes and ask user if selection is OK
%---------------------------------------------------------------------
axes(pH); hold on

while 1
    
    %------------------------------------------------------------------------------
    % Plot mesh
    %------------------------------------------------------------------------------
    patchinfo.Vertices          = [xnodes,ynodes];
    patchinfo.Faces             = t(ti,:);
    patchinfo.FaceColor         = [0.0060    0.4086    0.8828];
    patchinfo.FaceVertexCData   = ones(size(xnodes,1),1)*[0 0 1];
    patchinfo.EdgeColor         = 'k';
    patchinfo.linewidth         = 1.3;
    triPatch                    = patch(patchinfo);
    
    
    % gui.sb.setText('Continue?')
    
    % Ask the user if they like their selection
    msg = [...
        'Select Continue to proceed with the selection you made. '...
        'Select Re-Try to select a different region.'];
    choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'Continue','Re-Try','Quit'},'DefaultOption',1,'Icon','Warning');
    drawnow; pause(0.05);  % this innocent line prevents the Matlab hang
    
    switch choice
        
        case 'Continue'
            
            break
            
        case 'Re-Try'
            
            delete(triPatch)
            
            % Ask user to draw polygon around region
            bpoly = DrawBoundingPolygon(fig,pH,'r');
            
            % gui.sb.setText('Finding elements in selected region...')
            
            % Find verts in polygon
            verts = find(PointInPolygon(xnodes,ynodes,bpoly(:,1), bpoly(:,2)));
            
            % We already have triangulation representation
            
            % Find all elements attached to vertices
            ti = vertexAttachments(trep,verts); ti = [ti{:}]';
            
            % Remove elements containing vertices that do not fall within the selected
            % polygon.
            ti = ti(sum(ismember(t(ti,:),verts),2) == 3);
            
            % gui.sb.setText('Highlighting elements in selected region...')
            
        case 'Quit'
            
            delete(bH); % Delete full boundary
            delete(triPatch); % Delete patch
            
            % gui.sb.setText('Ready')            
            return
            
    end

end

%---------------------------------------------------------------------
% We want to constraint elements edges on the interior of full mesh
% determine what edges these are
%---------------------------------------------------------------------
% gui.sb.setText('Getting boundary information from sub domain...')
fd = dtFreeBoundary(MESH);    % Free boundary edges of current mesh
sd = dtFreeBoundary(MESH,ti); % Free boundary edges of sub mesh

% Connect interior node strings to constrain
[sb,intStr] = compileNodeStrings(sd,fd); %#ok<ASGLU>

% clear up some memory
clear sb fd sd verts

%---------------------------------------------------------------------
% Extract external boundary
%---------------------------------------------------------------------
% To reduce the chance of bad connections from the external boundary to the
% interior constraint, append the interior node string to the bounding
% polygon.
IN = PointInPolygon(bpoly(:,1), bpoly(:,2),PTS.Poly(1).x,PTS.Poly(1).y);

if all(~IN); errordlg('Bad Selection','Error'); return; end

% Create a column vector representing indices in x & y
ix = (1:length(bpoly(:,1)))';

% Set indices that do not below to zero
ix(IN == 1) = 0;

% Shift indices so we can extract one column vector of points
shiftN = abs(find(ix == 0,1,'last') - length(ix));

% Shift
ix = circshift(ix,shiftN);

% Remove zeros
ix(ix == 0) = [];

% Store relevant boundary (partial boundary)
bpolyx = bpoly(ix,1);
bpolyy = bpoly(ix,2);


% Number of nodes in main string
n = numel(intStr(1).str);

% Detmerine how to concatenate partial boundary with node string
intC = [... % Interior path constraint
    xnodes(intStr(1).str([1,n])),...
    ynodes(intStr(1).str([1,n]))];

% Find closest end point to the first point in our partial boundary p
closestPoint = knnsearch(intC, [bpolyx(end), bpolyy(end)]);

if closestPoint == 1 % Append Interior path constraint to end
    
    % New domain nodes
    bpolyx = [bpolyx; xnodes(intStr(1).str)];
    bpolyy = [bpolyy; ynodes(intStr(1).str)];

    % Ensure we have a unique list of points
    bpoly = unique([bpolyx,bpolyy],'rows','stable');
    
    % Close polygon
    bpoly(end+1,:) = bpoly(1,:);
    
else % Flip the path and append to end
    
    % New domain nodes
    bpolyx = [bpolyx; flipud(xnodes(intStr(1).str))];
    bpolyy = [bpolyy; flipud(ynodes(intStr(1).str))];

    % Ensure we have a unique list of points
    bpoly = unique([bpolyx,bpolyy],'rows','stable');
    
    % Close polygon
    bpoly(end+1,:) = bpoly(1,:);
    
end

% hold on
% h1 = plot(bpoly(:,1), bpoly(:,2),'-r*');
% 
% pause(5);
% delete(h1)

%------------------------------------------
% Repeat the above for boundary extraction
%------------------------------------------

% Use original bounding box to find points to use
[IN,~,~] = InPolygon(PTS.Poly(1).x,PTS.Poly(1).y,bpoly(:,1), bpoly(:,2));
%IN = PointInPolygon(PTS.Poly(1).x,PTS.Poly(1).y,bpoly(:,1), bpoly(:,2));

if all(~IN); errordlg('Bad Selection','Error'); return; end

% Create a column vector representing indices in x & y
ix = (1:length(PTS.Poly(1).x))';

% Set indices that do not below to zero
ix(IN == 0) = 0;

% Shift indices so we can extract one column vector of points
shiftN = abs(find(ix == 0,1,'last') - length(ix));

% Shift
ix = circshift(ix,shiftN);

% Remove zeros
ix(ix == 0) = [];

% Store relevant boundary (partial boundary)
pb.xy = [PTS.Poly(1).x(ix),PTS.Poly(1).y(ix)];

% Number of nodes in main string
n = numel(intStr(1).str);

% Detmerine how to concatenate partial boundary with node string
intC = [... % Interior path constraint
    xnodes(intStr(1).str([1,n])),...
    ynodes(intStr(1).str([1,n]))];

% Find closest end point to the first point in our partial boundary p
closestPoint = knnsearch(intC,pb(1).xy(end,:));

% clear up some memory
clear IN ix n intC 

%---------------------------------------------------------------------
% Create temporary PTS
%---------------------------------------------------------------------
if closestPoint == 1 % Append Interior path constraint to end
    
    % New domain nodes
    xt = [pb(1).xy(:,1); xnodes(intStr(1).str)];
    yt = [pb(1).xy(:,2); ynodes(intStr(1).str)];

    % Ensure we have a unique list of points
    xy = unique([xt,yt],'rows','stable');
    
    % Close polygon
    xy(end+1,:) = xy(1,:);
    
    % Store in temporary structure
    tPTS.Poly(1).x = xy(:,1);
    tPTS.Poly(1).y = xy(:,2);
    
else % Flip the path and append to end
    
    % New domain nodes
    xt = [pb(1).xy(:,1); flipud(xnodes(intStr(1).str))];
    yt = [pb(1).xy(:,2); flipud(ynodes(intStr(1).str))];

    % Ensure we have a unique list of points
    xy = unique([xt,yt],'rows','stable');
    
    % Close polygon
    xy(end+1,:) = xy(1,:);
    
    % Store in temporary structure
    tPTS.Poly(1).x = xy(:,1);
    tPTS.Poly(1).y = xy(:,2);
    
end

% Treat constraint as external barrier
tPTS.Constraints = [];
tPTS.Constraints(1).num = 3;
tPTS.Constraints(1).type  = 'External Barrier';
tPTS.Constraints(1).xy    = unique([xnodes(intStr(1).str),ynodes(intStr(1).str)],'rows','stable');
tPTS.Constraints(1).data  = [];

% clear up some memory
clear xt yt xy pb intStr

%---------------------------------------------------------------------
% Determine if additional polygons should be included in temporary domain
%---------------------------------------------------------------------
for k = 2:length(PTS.Poly)
    
    IN = PointInPolygon(PTS.Poly(k).x,PTS.Poly(k).y,tPTS.Poly(1).x, tPTS.Poly(1).y);
    
    if all(IN) % Add
        
        tPTS.Poly(end+1).x = PTS.Poly(k).x;
        tPTS.Poly(end).y   = PTS.Poly(k).y;
        
    end
    
end

%---------------------------------------------------------------------
% Plot subdomain boundary
%---------------------------------------------------------------------
delete(bH); % Delete full boundary
delete(triPatch); % Delete patch

%----------------------------------------------------------------------
% Plot boundary
%----------------------------------------------------------------------

gui.sb.ProgressBar.setVisible(true)
gui.sb.ProgressBar.setIndeterminate(true);
% gui.sb.setText('Removing elements...')

% Remove faces from current mesh
meshPatch   = findobj(pH,'Tag','Mesh');
faces       = get(meshPatch(end),'faces');
it          = (1:size(t,1))'; it(ti,:) = [];
faces       = faces(it,:);
set(meshPatch(end) ,'faces',faces)

% Plot polygon
for k = 1:length(tPTS.Poly)
    
    hold on
    bH = plot(tPTS.Poly(k).x,tPTS.Poly(k).y,'linewidth',2,'Color',[.7 .5 0]);
    set(bH ,'Tag','sub domain');
    
end

% Plot subdomain constraint
for k = 1:length(tPTS.Constraints)
    
    hold on
    p = tPTS.Constraints(k).xy;
    bH = plot(p(:,1),p(:,2),'linewidth',2,'Color','r');
    set(bH ,'Tag','sub domain');
    
end


%---------------------------------------------------------------------
% Set admesh button to run sub meshing routine
%---------------------------------------------------------------------
% gui.elementsToRemove = ti;
% gui.subPTS = tPTS;
% % guidata(fig,gui);
% set(gui.Window , 'windowkeypressfcn' , @CancelADMESHSubMesh)
% set(gui.RunAdmeshButton, 'Callback'  , {@ADmeshSubMeshRoutine,fig})

uiwait(msgbox(['ADMESH is ready to mesh subdomain. Enter in the appropriate ',...
    'parameters and select ''Run ADMESH''. ADMESH treats your input ',...
    'differenly now. To return back to the normal settings, push the ''c'' ',...
    'key on you keyboard.'],'ADMESH','modal'));

% gui.sb.setText('Ready.')

end