classdef MESHobj
    
    % Create instance of a MESHobj
    %
    % Syntax
    %
    %       MESH = MESHobj([])
    %       MESH = MESHobj('filename')
    %       MESH = MESHobj(t,x,y)
    %       MESH = MESHobj(t,x,y,z)
    %       MESH = MESHobj(t,x,y,z,C)
    %
    % Description
    %
    %     MESHobj creates an instance of the mesh class, which contains
    %     information on the mesh.
    %
    % Author: Dustin West
    % The Ohio State University
    % email address: dww.425@gmail.com
    % October 2013; Last revision: 21-October-2013
    
    properties
        
        %Points: The coordinates of the points in the triangulation
        Points
        
        %ConnectivityList: The computed triangulation connectivity list
        ConnectivityList
        
        %Constraints; The imposed edge constraints data structure
        Constraints
        
    end
    
    methods
        
        function MESH = MESHobj(varargin)
            
            
            % 1 Input argument
            if nargin == 1
                
                % Check for input type
                if isa(varargin{1},'double')
                    
                    % Return empty mesh object
                    return
                    
                else
                    
                    % Get file name
                    filename = varargin{1};
                    
                    [t,p,c] = Read14File(file);
                    
                end
                
                
                
                
            end
            
            
        end
    end
    
end

function [t,p,c] = Read14File(file)
% Mesh2EdgeStructure - Reads in domain and it's attributes
%
% Syntax:  [mesh] = Read14File(file)
%
% Inputs:
%    guiFig - GUI Object Handle
%    guiH   - GUI Data Structure
%
% Outputs:
%    mesh
%
% See also: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% October 2013; Last revision: 21-October-2013

%------------------------------- BEGIN CODE -----------------------------------

%------------------------------------------------------------------------------
% Open grid file and read in nodes, coordinates and connectivity.
%------------------------------------------------------------------------------
fid=fopen(file, 'r');    % Open file

fgetl(fid);              % Read in Grid name

uiStatusBar('Reading in nodes & coordinates...')

% Read in the number of elements & number of grid points
C = textscan(fid,'%f %f %f %*[^\n]',1);

% Check for element connectivity list
if isnan(C{3}) || C{3} == 3
    
    Nn = 3;
    
else
    
    Nn = 6;
    
end

%------------------------------------------------------------------------------
% Initialize MESH Structure
%------------------------------------------------------------------------------
MESH = struct('Points',zeros(C{2},3),'ConnectivityList',zeros(C{1},3),'cpplon',[],'cpplat',[]);

%------------------------------------------------------------------------------
% Read in node numbers and coordinates
%------------------------------------------------------------------------------

% If the file is large, read in blocks for speed up
blocks = 10000; Nc = C{2};

if Nc > blocks
    
    p = cell(ceil(Nc/blocks),1); n = 1;
    
    while 1
        
        p{n} = cell2mat(textscan(fid, '%*f %f %f %f %*[^\n]', blocks));
        
        Nc = Nc - blocks;
        
        n = n+1;
        
        if Nc < blocks
            p{n} = cell2mat(textscan(fid, '%*f %f %f %f %*[^\n]', Nc));
            break
        end
        
    end
    
    p = cell2mat(p);
    
else
    
    p = cell2mat(textscan(fid, '%*f %f %f %f %*[^\n]', C{2}));
    
end

% Assign point list to mesh
MESH.Points = p; clear p

% Reverse sign on elevation in mesh
MESH.Points(:,3) = -MESH.Points(:,3);

%------------------------------------------------------------------------------
% Read in the element connectivity table
%------------------------------------------------------------------------------
uiStatusBar('Reading in element connectivity table...')

% If the file is large, read in blocks for speed up
blocks = 10000; Ne = C{1};

if Ne > blocks
    
    conn = cell(ceil(Ne/blocks),1); n = 1;
    
    while 1
        
        if Nn == 3
            conn{n} = cell2mat(textscan(fid, '%*f %*f %f %f %f %*[^\n]', blocks));
        else
            conn{n} = cell2mat(textscan(fid, '%*f %*f %f %f %f %f %f %f %*[^\n]', blocks));
        end
        
        Ne = Ne - blocks;
        
        n = n+1;
        
        if Ne < blocks
            if Nn == 3
                conn{n} = cell2mat(textscan(fid, '%*f %*f %f %f %f %*[^\n]', Ne));
            else
                conn{n} = cell2mat(textscan(fid, '%*f %*f %f %f %f %f %f %f %*[^\n]', Ne));
            end
            break
        end
        
    end
    
    conn = cell2mat(conn);
    
else
    
    if Nn == 3
        conn = cell2mat(textscan(fid, '%*f %*f %f %f %f %*[^\n]', Ne));
    else
        conn = cell2mat(textscan(fid, '%*f %*f %f %f %f %f %f %f %*[^\n]', Ne));
    end
    
    
end

if Nn == 3
    
    % Assign connectivity to mesh
    MESH.ConnectivityList = conn; clear conn
    
    % Create element connectivity
    trep = triangulation(MESH.ConnectivityList,MESH.Points(:,[1,2]));
    
    MESH.ElementConnectivity = neighbors(trep);
    
    clear trep
    
else
    
    % Assign connectivity to mesh
    MESH.ConnectivityList = conn(:,[1 2 3]);
    
    % Assign element connectivity to mesh
    MESH.ElementConnectivity = conn(:,[4 5 6]);
    
    clear conn
    
end


%------------------------------------------------------------------------------
% Get coordinate information for file. Convert coordinates, if needed
%------------------------------------------------------------------------------

% Determine if coordinates are in geographic range
if all(MESH.Points(:,1) < 180 & MESH.Points(:,1) > -180 & MESH.Points(:,2) < 90 & MESH.Points(:,2) > -90)
    
    % Construct a questdlg with three options
    msg = 'Select the coordinate system the mesh file is in.';
    choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'Geographic','Meters','None'},'DefaultOption',1,'Icon','Warning');
    drawnow; pause(.0005);
    
    MESH.cpplon = [];
    MESH.cpplat = [];
    
    % Handle response
    switch choice
        
        case 'Geographic'
            
            uiStatusBar('Converting to cartesian coordinate system...')
            
            % Convert to XY (meters)
            [MESH.Points(:,1),MESH.Points(:,2),MESH.cpplon,MESH.cpplat] = ...
                Geo2Meters(MESH.Points(:,1),MESH.Points(:,2));
            
    end
    
end

%------------------------------------------------------------------------------
% Read in open ocean boundary segments.
%------------------------------------------------------------------------------

% Read in the number of open ocean segments
NOPE = cell2mat(textscan(fid, '%f %*[^\n]', 1)); if isempty(NOPE); NOPE = 0; end

% Scan over the total number of elevation specified boundary nodes
textscan(fid, '%*f %*[^\n]', 1);

% Store open ocean boundary coordinates
if NOPE ~= 0
    
    uiStatusBar('Reading in open ocean boundary segements...')
    
    % Intitialize PTS Constraint structure
    MESH.Constraints(NOPE,1) = struct('num',[],'type',[],'nodeStr',[],'data',[]);
    
    for k = 1:NOPE
        
        % Tag for Open Ocean Boundary
        MESH.Constraints(k).num = -1;
        
        MESH.Constraints(k).type = 'Open';
        
        % Read number of nodes in open ocean segment k
        NVDLL = cell2mat(textscan(fid, '%f %*[^\n]', 1));
        
        % Node numbers in boundary segment k
        Nodes = cell2mat(textscan(fid, '%f %*[^\n]', NVDLL));
        
        % Store node string of Open Ocean Boundary
        MESH.Constraints(k).nodeStr = Nodes;
        
    end
    
end

%------------------------------------------------------------------------------
% Read in normal flow specified boundary segements.
%------------------------------------------------------------------------------

% Read in the number of normal flow specified boundary segments
NBOU = cell2mat(textscan(fid, '%f %*[^\n]', 1)); if isempty(NBOU); NBOU = 0; end

% Scan over the total number of normal flow specified boundary nodes
textscan(fid, '%*f %*[^\n]', 1);

% Store normal flow boundary coordinates
if NBOU ~= 0
    
    uiStatusBar('Reading in normal flow specified boundary segements...')
    
    % Add space to PTS data structure
    if isfield(MESH,'Constraints')
        MESH.Constraints(end+NBOU) = struct('num',[],'type',[],'nodeStr',[],'data',[]);
    else
        MESH.Constraints(NBOU,1)   = struct('num',[],'type',[],'nodeStr',[],'data',[]);
    end
    
    for k = (NOPE+1):(NOPE+NBOU) % Loop over all the normal flow specified boundary segments
        
        C = textscan(fid, '%f %f %*[^\n]', 1); [NVELL, BoundType] = deal(C{:});
        
        if sum(BoundType == [0 2 10 12 20 22 30]) % Read in external boundary data
            
            % Tag External Boundary
            MESH.Constraints(k).type = 'External Boundary';
            
            % Store value
            MESH.Constraints(k).num = BoundType;
            
            % Node numbers in boundary segment k
            Nodes = cell2mat(textscan(fid, '%f %*[^\n]', NVELL));
            
            % Node String
            MESH.Constraints(k).nodeStr = Nodes;
            
        elseif sum(BoundType == [1 11 21]) % Read in internal boundary data
            
            % Tag Internal Boundary
            MESH.Constraints(k).type = 'Internal Boundary';
            
            % Store value
            MESH.Constraints(k).num = BoundType;
            
            % Node numbers in boundary segment k
            Nodes = cell2mat(textscan(fid, '%f %*[^\n]', NVELL));
            
            % Node String
            MESH.Constraints(k).nodeStr = Nodes;
            
        elseif sum(BoundType == [3 13 23])
            
            % Tag External Barrier
            MESH.Constraints(k).type = 'External Barrier';
            
            % Store value
            MESH.Constraints(k).num = BoundType;
            
            % Node numbers in boundary segment k
            Nodes = cell2mat(textscan(fid, '%f %f %f %*[^\n]', NVELL));
            
            % Node String
            MESH.Constraints(k).nodeStr = Nodes(:,1);
            
            % Data - [barrier height, coefficient]
            MESH.Constraints(k).data = Nodes(:,2:end);
            
        elseif sum(BoundType == 18)
            
            % Tag Channel Line
            MESH.Constraints(k).type = 'Line';
            
            % Store value
            MESH.Constraints(k).num = BoundType;
            
            % Node numbers in boundary segment k
            Nodes = cell2mat(textscan(fid, '%f %f %*[^\n]', NVELL));
            
            % Node String
            MESH.Constraints(k).nodeStr = Nodes(:,1);
            
            % Data - [barrier height, coefficient]
            MESH.Constraints(k).data = Nodes(:,2);
            
        elseif sum(BoundType == [4 24]) % Read in internal barrier boundary data
            
            % Tag Internal Barrier
            MESH.Constraints(k).type = 'Internal Barrier Type 1';
            
            % Store value
            MESH.Constraints(k).num = BoundType;
            
            % Node numbers in boundary segment k
            Nodes = cell2mat(textscan(fid, '%f %f %f %f %f %*[^\n]', NVELL));
            
            % Node String
            MESH.Constraints(k).nodeStr = Nodes(:,1:2);
            
            % Data - [barrier height, coefficient]
            MESH.Constraints(k).data = Nodes(:,3:end);
            
        elseif sum(BoundType == [5 25]) % Read in internal barrier boundary data
            
            % Tag Internal Barrier
            MESH.Constraints(k).type = 'Internal Barrier Type 2';
            
            % Store value
            MESH.Constraints(k).num = BoundType;
            
            % Node numbers in boundary segment k
            Nodes = cell2mat(textscan(fid, '%f %f %f %f %f %f %f %f %*[^\n]', NVELL));
            
            % Node String
            MESH.Constraints(k).nodeStr = Nodes(:,1:2);
            
            % Data - [barrier height, coefficient]
            MESH.Constraints(k).data = Nodes(:,3:end);
            
        end
        
    end
    
end

if ~isfield(MESH,'Constraints')
    MESH.Constraints = [];
end

fclose(fid); % Close the file

uiStatusBar('File read in complete')

end