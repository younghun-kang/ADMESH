function [MESH,xyzFun,status] = Read14File(file,app)

%% READ14FILE Reads in a fort.14 mesh file
%   Fort14 = Read14File() reads in a user-selected fort.14 file and saves
%   it as a structure with the following fields:
%
%     Fort14.MeshID: the mesh identification string
%     Fort14.TR: the mesh saved as a MATLAB?triangulation object
%     Fort14.Boundary: mesh boundary saved as a MATLAB?polyshape object. 
%     Fort14.BoundarySegments: a 1 by N structure, where N is the number
%       of boundary segments, with fields:
%          BoundarySegments(i).Type = string identification of the
%            boundary segment type, either 'Elevation' or 'NormalFluxTypeK'
%            where K is an integer value defining the particular normal-
%            flux boundary type; see http://adcirc.org/ for details.
%          BoundarySegments(i).Nodes = node numbers that make up boundary
%            segment i.
%          BoundarySegments(i).Data = additional node-based boundary data
%            associated with some normal-flux boundary types; see
%            http://adcirc.org/ for details. If no additional data is
%            associated with a boundary type, the field is empty.
%
%   Fort14 = Read14File(FileName) same as above but reads in the fort.14
%   file specified by FileName.
                                
%
%   NOTE: See the ADCIRC documention at https://adcirc.org/ for more
%   details on the format of the input/output files.
%
%   See also Read63File, Read64File, Plot14File
                    
                          
%

%% Open the fort.14 file ==================================================
try
status = 1;
fileID = fopen(file,'r');

%% Create the Fort14 structure ============================================

MESH = struct('MeshID',[],'TR',[],'Boundary',[],'BoundarySegments',[]);

%% Read in the mesh =======================================================

% Read in the header info
MESH.MeshID = fgetl(fileID);

% Read in the number of elements and the number of nodes
NumberOfElems = fscanf(fileID,'%u',1);
NumberOfNodes = fscanf(fileID,'%u %* [^\n]',1);

% Read in the point list
Points = textscan(fileID,'%f %f %f %f %*[^\n]',NumberOfNodes);

% Read in the element connectivity
ConnectivityList = cell2mat(textscan(fileID,'%f %f %f %f %f %*[^\n]',...
    NumberOfElems));

% Contruct the Matlab triangulation
MESH.TR = triangulation(ConnectivityList(:,3:5),Points{2},Points{3},...
    Points{4});

%% Read in the boundary data ==============================================

% Read in the number of elevation-specified boundary segments
ElevSegments = cell2mat(textscan(fileID,'%u %*[^\n]',2));

% Loop over the elevation-specified boundary segments
ExternalBoundarySegments = {};
if ~isempty(ElevSegments)
    for i = 1:ElevSegments(1)                                         
        Segment = textscan(fileID,'%f %*[^\n]',1);
        MESH.BoundarySegments(i).Type = 'Elevation';
        MESH.BoundarySegments(i).Nodes = fscanf(fileID,'%u',Segment{1});
        MESH.BoundarySegments(i).Data = [];
        ExternalBoundarySegments{i} = MESH.BoundarySegments(i).Nodes;
    end
else
    ElevSegments = 0;
end

% Read in the number of flux-specified boundary segments
FluxSegments = cell2mat(textscan(fileID,'%u %*[^\n]',2));

% Loop over the flux-specified boundary segments
X = {};
Y = {};
if ~isempty(FluxSegments)
    j = ElevSegments(1);

    k = 0;
    for i = 1:FluxSegments(1)
        Segment = textscan(fileID,'%f %f %*[^\n]',1);
        MESH.Constraints(i).num = Segment{2};
        MESH.BoundarySegments(i).Type = ...
            ['NormalFluxType',num2str(Segment{2})];
        switch Segment{2}
            case {0, 1, 2, 10, 11, 12, 20, 21, 22, 30} % No-normal flow ===
                SegmentData = textscan(fileID,'%u %*[^\n]',Segment{1});
                MESH.BoundarySegments(i).Nodes = SegmentData{1};
                MESH.Constraints(i).nodeStr(:,1) = SegmentData{1};
                MESH.BoundarySegments(i).Data = [];               
                switch Segment{2}                    
                    case {0, 2, 10, 12, 20, 22, 30} % External ============
                        % If the boundary segment is actually an internal
                        % boundary,....
                        if MESH.BoundarySegments(i).Nodes(1) == ...
                                MESH.BoundarySegments(i).Nodes(end) ...
                                &&  ElevSegments(1)~=0                                
                            % Then change the type and save it as a "hole"
                            % in the boundary polygon
                            MESH.BoundarySegments(i).Type = ...
                                ['NormalFluxType',num2str(Segment{2}+1)];
                            k = k + 1;
                            X{k} = MESH.TR.Points(unique(MESH.BoundarySegments(i).Nodes,'stable'),1);
                            Y{k} = MESH.TR.Points(unique(MESH.BoundarySegments(i).Nodes,'stable'),2);
                        else
                            j = j + 1;
                            ExternalBoundarySegments{j} = MESH.BoundarySegments(i).Nodes;                
                        end
                    case {1, 11, 21 } % Internal ==========================
                        k = k + 1;
                        X{k} = MESH.TR.Points(unique(MESH.BoundarySegments(i).Nodes,'stable'),1);
                        Y{k} = MESH.TR.Points(unique(MESH.BoundarySegments(i).Nodes,'stable'),2);
                                                                     
                end
            case {3, 13, 23} % External barriers ==========================                
                SegmentData = textscan(fileID,'%u %f %f %*[^\n]',Segment{1});
                MESH.BoundarySegments(i).Nodes = SegmentData{1};
                MESH.BoundarySegments(i).Data = [SegmentData{2},SegmentData{3}];
                MESH.Constraints(i).nodeStr(:,1) = SegmentData{1};
                j = j + 1;
                ExternalBoundarySegments{j} = MESH.BoundarySegments(i).Nodes;                                            
            case {18,6,19,17}
                SegmentData = textscan(fileID,'%u %f %f %f %*[^\n]',Segment{1});
                MESH.BoundarySegments(i).Nodes = SegmentData{1};
                MESH.Constraints(i).nodeStr(:,1) = SegmentData{1};
                MESH.BoundarySegments(i).Data = [SegmentData{2},SegmentData{3},...
                    SegmentData{4}];
                MESH.Constraints(i).data = [SegmentData{2},SegmentData{3},...
                    SegmentData{4}];
            case {4, 24} % Internal barriers ==============================               
                SegmentData = textscan(fileID,'%u %f %f %f %f %*[^\n]',Segment{1});
                MESH.BoundarySegments(i).Nodes = SegmentData{1};
                MESH.Constraints(i).nodeStr(:,1) = SegmentData{1};
                MESH.BoundarySegments(i).Data = [SegmentData{2},SegmentData{3},...
                    SegmentData{4},SegmentData{5}];
                k = k + 1;
                X{k} = [ MESH.TR.Points(MESH.BoundarySegments(i).Nodes,1);...
                   flipud(MESH.TR.Points(MESH.BoundarySegments(i).Data(:,1),1)) ];    
                Y{k} = [ MESH.TR.Points(MESH.BoundarySegments(i).Nodes,2);...
                   flipud(MESH.TR.Points(MESH.BoundarySegments(i).Data(:,1),2)) ]; 
               % Note: internal boundary segments may intersect with the
               % outer land boundary in which case they should be saved as
               % external boundary segments for purposes of constructing
               % the polyshape of the boundary (see code below)
               j = j + 1;
               ExternalBoundarySegments{j} = ...
                   [ MESH.BoundarySegments(i).Nodes(1),...
                     MESH.BoundarySegments(i).Data(1,1) ];
               j = j + 1;  
               ExternalBoundarySegments{j} = ...
                   [ MESH.BoundarySegments(i).Nodes(end),...
                     MESH.BoundarySegments(i).Data(end,1) ];   
            case {5, 25} % Internal barrier with pipe =====================
                SegmentData = textscan(fileID,'%u %f %f %f %f %f %f %f %*[^\n]',Segment{1});
                MESH.Constraints(i).nodeStr(:,1) = SegmentData{1};
                MESH.BoundarySegments(i).Nodes = SegmentData{1};
                MESH.BoundarySegments(i).Data = [SegmentData{2},SegmentData{3},...
                    SegmentData{4},SegmentData{5},SegmentData{6},SegmentData{7},SegmentData{8}];
                k = k + 1;
                X{k} = [ MESH.TR.Points(MESH.BoundarySegments(i).Nodes,1);...
                   flipud(MESH.TR.Points(MESH.BoundarySegments(i).Data(:,1),1)) ];
                Y{k} = [ MESH.TR.Points(MESH.BoundarySegments(i).Nodes,2);...
                   flipud(MESH.TR.Points(MESH.BoundarySegments(i).Data(:,1),2)) ];
        end
    end
end
%% Construct a polyshape of the boundary ==================================
if length(ExternalBoundarySegments) > 0
ExternalBoundaryNodes = ExternalBoundarySegments{1};
ExternalBoundarySegments{1} = [];
First = ExternalBoundaryNodes(1);
Last  = ExternalBoundaryNodes(end);
while First ~= Last                                                     
    FLAG = 0;
    for i = 1:j
        if ~isempty(ExternalBoundarySegments{i})
            if Last == ExternalBoundarySegments{i}(1)
                ExternalBoundaryNodes = [ ExternalBoundaryNodes; ...
                    ExternalBoundarySegments{i}(2:end) ];
                ExternalBoundarySegments{i} = [];
                FLAG = 1;
            elseif Last == ExternalBoundarySegments{i}(end)
                ExternalBoundaryNodes = [ ExternalBoundaryNodes; ...
                    flipud(ExternalBoundarySegments{i}(1:end-1)) ];
                ExternalBoundarySegments{i} = [];
                FLAG = 1;
            elseif First == ExternalBoundarySegments{i}(1)
                ExternalBoundaryNodes = [
                    flipud(ExternalBoundarySegments{i}(2:end));...
                    ExternalBoundaryNodes ];
                ExternalBoundarySegments{i} = [];
                FLAG = 1;
            elseif First == ExternalBoundarySegments{i}(end)
                ExternalBoundaryNodes = [ 
                    ExternalBoundarySegments{i}(1:end-1);...
                    ExternalBoundaryNodes ];
                ExternalBoundarySegments{i} = [];
                FLAG = 1;
            end
        end
        First = ExternalBoundaryNodes(1);
        Last = ExternalBoundaryNodes(end);
    end
    
    if FLAG == 0
        ExternalBoundaryNodes = [];
        break;
    end
end
ExternalBoundaryNodes = ExternalBoundaryNodes(1:end-1);
OuterX = MESH.TR.Points(ExternalBoundaryNodes,1);
OuterY = MESH.TR.Points(ExternalBoundaryNodes,2);

% Construct outer boundary polyshape (oriented ccw)
warning('off','all');
OuterBoundary = polyshape(OuterX,OuterY,'SolidBoundaryOrientation','ccw');
% Construct inner boundary polyshape (oriented cw)
if ~isempty(X)
InnerBoundary = polyshape(X,Y,'SolidBoundaryOrientation','cw');
% Combine outer and inner boundaries to make one polyshape
MESH.Boundary = addboundary(OuterBoundary,InnerBoundary.Vertices);
end
end
% Close file
fclose(fileID);

MESH.Points = MESH.TR.Points;
MESH.ConnectivityList = MESH.TR.ConnectivityList;

if all(abs(MESH.Points(:,3)) == 1) || all(MESH.Points(:,3) == 0)
        
        xyzFun = [];
        
    else
        % Turn on progress bar
        %         gui.sb.ProgressBar.setVisible(true)
        %
        %         set(gui.sb.Progressbar,'Minimum',0,'Maximum',5,'Value',1)
        %
        %         xmin = min(MESH.Points(:,1));
        %         xmax = max(MESH.Points(:,1));
        %         ymin = min(MESH.Points(:,2));
        %         ymax = max(MESH.Points(:,2));
        %
        %         Offsetx = (xmax-xmin)*.05;
        %         Offsety = (ymax-ymin)*.05;
        %
        %         % Generate unique x and y grid vectors
        %         xRange = unique(MESH.Points(:,1)');
        %         yRange = unique(MESH.Points(:,2)');
        %
        %         % Produce a grid that won't run out of memory
        %         maxVecLength = 1000;
        %
        %         % Space vectors accordingly
        %         if numel(xRange) > maxVecLength
        %             spacing = floor(numel(xRange)/maxVecLength);
        %             xRange = xRange(1:spacing:end);
        %         end
        %
        %         if numel(yRange) > maxVecLength
        %             spacing = floor(numel(yRange)/maxVecLength);
        %             yRange = yRange(1:spacing:end);
        %         end
        %
        %         xRange = [xmin-Offsetx, xRange, xmax+Offsetx];
        %         yRange = [ymin-Offsety, yRange, ymax+Offsety];
        %
        %         % Create grid
        %         [xq,yq] = meshgrid(xRange,yRange);
        %
        %         set(gui.sb.Progressbar,'Minimum',0,'Maximum',5,'Value',2)
        %
        %         % Interpolate scatter set
        %         zq = griddata(MESH.Points(:,1),MESH.Points(:,2),MESH.Points(:,3),xq,yq);
        %
        %         set(gui.sb.Progressbar,'Minimum',0,'Maximum',5,'Value',4)
        %
        %         % Interpolate NaN values via inpaint method
        %         zq = inpaint_nans(zq,4);
        %
        %         % Create gridded interpolant
        %         xyzFun = griddedInterpolant(xq',yq',zq','linear','nearest');
        %
        %         set(gui.sb.Progressbar,'Minimum',0,'Maximum',5,'Value',5)
        %
        %         % Turn off progress bar
        %         gui.sb.ProgressBar.setVisible(false)
                        
        msg = 'Creating elevation interpolant function...';
        progdlg = uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');
        
        xyzFun = scatteredInterpolant(MESH.Points(:,1),MESH.Points(:,2),MESH.Points(:,3),'linear','nearest');
        
        close(progdlg);
        
end
    
    
catch
    
    MESH   = [];    % MESH Structure
    xyzFun = [];    % Elevation
    status = 0;     % Error flag
    
end

