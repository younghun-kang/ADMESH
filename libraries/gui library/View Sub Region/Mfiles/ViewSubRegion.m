function ViewSubRegion(varargin)
% View Sub Region - Keeps specfied mesh in domain
%
% Syntax: ViewSubRegion(varargin)
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
fig = varargin{1}; gui = % guidata(fig);

% Assign variables
PTS     = gui.PTS;
MESH    = gui.MESH;

%---------------------------------------------------------------------
% Check if any data is loaded
%---------------------------------------------------------------------
pH = gui.ViewAxes; % Plot Window handle

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

% Get view status
status = get(varargin{1},'checked');

switch status
    
    case 'off' % Turn on
        
        set(varargin{1},'checked','on')
        
        
    case 'on'
        
        set(varargin{1},'checked','off')
        
        return
        % Return original mesh
        
        % Original mesh
        MESH = gui.originalMESH;
        
        %--------------------------------------------------------------------
        %Replace ADMESH MESH variable
        %--------------------------------------------------------------------
        %gui.MESH.ConnectivityList    = MESH.ConnectivityList(ti,:);
        %gui.subConnectivityListIndex = ti;
        
        
        
end

% Catch no PTS
% if ~isempty(MESH) && isempty(PTS)
%     
%     if ~isempty(MESH.Constraints)
%         PTS = Constraints2PTS(MESH);
%     else
%         PTS = triangulation2PTS(MESH);
%     end
%     
%     gui.PTS = PTS;
%     % guidata(gui.Window,gui)
%     
% end

%---------------------------------------------------------------------
% Select mesh domain
%---------------------------------------------------------------------

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

%---------------------------------------------------------------------
% Highlight nodes and ask user if selection is OK
%---------------------------------------------------------------------
axes(pH); hold on

while 1
    
    %-----------------------------------------------------------------
    % Plot mesh
    %-----------------------------------------------------------------
    patchinfo.Vertices          = [xnodes,ynodes];
    patchinfo.Faces             = t(ti,:);
    patchinfo.FaceColor         = [0.0343    0.5966    0.8199];
    patchinfo.EdgeColor         = 'k';
    patchinfo.linewidth         = 1.3;
    triPatch                    = patch(patchinfo);
    
    
    % gui.sb.setText('Continue?')
    
    % Ask the user if they like their selection
    msg = [
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
            
            % Find verts in polygon
            verts = find(PointInPolygon(xnodes,ynodes,bpoly(:,1), bpoly(:,2)));
            
            % We already have triangulation representation
            
            % Find all elements attached to vertices
            ti = vertexAttachments(trep,verts); ti = [ti{:}]';
            
            % Remove elements containing vertices that do not fall within the selected
            % polygon.
            ti = ti(sum(ismember(t(ti,:),verts),2) == 3);
            
            
        case 'Quit'
            
             set(varargin{1},'checked','off')
            delete(triPatch); % Delete patch
            % gui.sb.setText('Ready')      
            return
            
    end

end

%--------------------------------------------------------------------
%Replace ADMESH MESH variable
%--------------------------------------------------------------------

% Create unique node list
nodeList = unique(MESH.ConnectivityList(ti,:));

% Keep relevant mesh points
gui.MESH.Points = MESH.Points(nodeList,:);

% Keep relevant connectivity list
gui.MESH.ConnectivityList = MESH.ConnectivityList(ti,:);
[~,gui.MESH.ConnectivityList] = ismember(gui.MESH.ConnectivityList,nodeList);

gui.subConnectivityListIndex = ti;

gui.MESH.Constraints = [];

%--------------------------------------------------------------------
% Store original mesh
%--------------------------------------------------------------------

% Remove section from connectivity List
MESH.ConnectivityList(ti,:) = [];

% Store mesh
gui.originalMESH = MESH;

PlotMesh(gui,.05)

% guidata(gui.Window,gui)

% gui.sb.setText('Ready.')


end