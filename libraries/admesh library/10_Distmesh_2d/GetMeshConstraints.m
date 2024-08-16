function [p,nC,C,mesh] = GetMeshConstraints(p,hmin,PTS)
% ApplyMeshConstraints - Initializes variables for constrained
% triangulations
%
% Syntax:  [p,nC,C,mesh] = ApplyMeshConstraints(hmin,guiFig)
%
% Inputs:
%    p      - point list
%    hmin   - Minumum elements size
%    guiFig - handle that identifies the figure
%
% Outputs:
%    mesh  - structure containing constraints
%    pfix  - coordinate list of fixed nodes
%    nfix  - number of coordinate pairs in pfix
%    C     - list of edge constraints that are defined by
%            an numc-by-2 matrix, numc being the number of constrained
%            edges
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% August 2013; Last revision: 08-August-2013

%--------------------------------------------------------------------------
% Check for constraints
%--------------------------------------------------------------------------
if isempty(PTS.Constraints) % No constraints exist
    
    % Set mesh struct to []
    mesh.Constraints = [];
    
    % Set delaunay constraints matrix to empty
    C = [];
    
    % 0 points to constrain
    nC = 0;
    
    % remove duplicate nodes
    p=unique(p,'rows','stable');
    
    return
    
end

%--------------------------------------------------------------------------
% Determine the total number of constraint segments (excluding ocean BC)
%--------------------------------------------------------------------------
% Number designation for constraints
cNum = [3 13 23 4 5 24 25 18 19 17];

% Location of constraints in PTS
idc = find(ismember([PTS.Constraints.num],cNum))';

% Number of segments
ncSeg = numel(idc);

%--------------------------------------------------------------------------
% Determine the total number of ocean BC segments 
%--------------------------------------------------------------------------
% Number designation for ocean BC
cNum = -1;

% Location of ocean BC in PTS
ido = find(ismember([PTS.Constraints.num],cNum))';

% Number of ocean BC
noSeg = numel(ido);

%--------------------------------------------------------------------------
% Double check incase constraints field exist but with no constraints
%--------------------------------------------------------------------------
if isempty(idc) && isempty(ido) 
    
    % Set mesh struct to []
    mesh.Constraints = [];
    
    % Set delaunay constraints matrix to empty
    C = [];
    
    % 0 points to constrain
    nC = 0;
    
    % remove duplicate nodes
    p=unique(p,'rows','stable');
    
    return
    
end

%--------------------------------------------------------------------------
% Create Constraints Matrix
%--------------------------------------------------------------------------

% Initialize cell, pC
pC = cell(ncSeg + noSeg,1);

% Initialize mesh structure
mesh.Constraints(ncSeg + noSeg,1) = ...
    struct('num',[],'type',[],'nodeStr',[],'data',[]);

% Counters
i = 1; j = 1;

% End index in constraints
e   = 0;

if isempty(idc) % No constraint segments, only ocean BC
        
    % Set delaunay constraints matrix to empty
    C = [];
    
else
    
    % Initialize delaunay constraints matrix as a cell
    C = cell(ncSeg,1);
           
    % Build constraints matrix by segmenting each constraint.
    
    %--------------------------------------------------------
    % Line Constraints
    %--------------------------------------------------------
    
    % delaunay triangulation requires a unique set of coordinates.
    % This poses a challenge for handling constraints with junctions
    % that have a common point. The following algorithm creates a
    % unique set of coordinates and modifies the constraints matrix
    % to handle junctions appropriately
    
    % Find all line constraints in PTS
    idx = find(ismember([PTS.Constraints.num],[19 18 17]))';
    
    if ~isempty(idx)
        
        % Create cell array for storing length of each line segment
        l = cell(numel(idx),1);
        
        % Loop over each line segment
        for k = 1:numel(idx)
            
            % Store coordinates in cell array. Space array equally based on
            % minimum element size
%             pC{i} = SpacePolyEqually(PTS.Constraints(idx(k)).xy,hmin);

            %--------------------------------------------------------------
            % 2021-03-18 Younghun: Enforce given points instead of uniform distribution.
            %--------------------------------------------------------------
            pC{i} = PTS.Constraints(idx(k)).xy;
            
%             %--------------------------------------------------------------
%             % 2021-04-12 Younghun: Empty points (the points will be set later using magnetic force.
%             % Uncomment this block for using magnetic force approach.
%             %--------------------------------------------------------------
%             pC{i} = [];
            
            % Store length of segment
            l{k} = size(pC{i},1);
            
            i = i + 1;
            
        end
        
        % Convert cell to matrix. Compute cumulative sum of the lengths
        pl = cell2mat(pC); l = cumsum(cell2mat(l));
        
        % Create unique coordinate set with back pointers
        [~,~,ic] = unique(pl,'rows','stable');
        
        % Determine the constraints matrix based on unique data set
        for k = 1:numel(idx)
            
            % Increase index
            s = e + 1; 
            e = l(k);
            
            C{j} = [ic(s:e-1) ic(s+1:e)];
            
            % Assign node string in mesh structure
            mesh.Constraints(j).nodeStr = ic(s:e);
            
            % Assign number associated with constraint
            mesh.Constraints(j).num = PTS.Constraints(idx(k)).num;
            
            % Assign node string type
            mesh.Constraints(j).type = PTS.Constraints(idx(k)).type;
            
            % Check for channel width data
            if ~isempty(PTS.Constraints(idx(k)).data)
                
                if isequal(PTS.Constraints(idx(k)).xy,pC{k})
                    mesh.Constraints(j).data = PTS.Constraints(idx(k)).data;
                else
                for jj = 1 : size(PTS.Constraints(idx(k)).data,2)
                % Create xyz scatter set
                xyz = [PTS.Constraints(idx(k)).xy PTS.Constraints(idx(k)).data(:,jj)];
                
                % Interpolate channel width to new coordinates
                F = scatteredInterpolant(xyz(:,1),xyz(:,2),xyz(:,3));
                
                mesh.Constraints(j).data(:,jj) = F(pC{k}(:,1),pC{k}(:,2));
                end
                end
                
            else
                
                mesh.Constraints(j).data  = zeros(size(pC{k},1),1);
                
            end
            
            j = j + 1;
            
            drawnow; % for graphics
            
        end
        
    end
        
    %--------------------------------------------------------
    % External barrier
    %--------------------------------------------------------
    
    % Find all external barrier constraints in PTS
    idx = find(ismember([PTS.Constraints.num],[3 13 23]))';
    
    if ~isempty(idx)

        for k = 1:numel(idx)
                        
            pC{i} = PTS.Constraints(idx(k)).xy;  % Assign pfix
            
            mesh.Constraints(j).data = PTS.Constraints(idx(k)).data;
            
            mesh.Constraints(j).num  = PTS.Constraints(idx(k)).num;
            
            mesh.Constraints(j).type = PTS.Constraints(idx(k)).type;
            
            % Increase index
            s = e + 1;
            e = size(pC{i},1) + e;

            % Do not want to constrain last point with first point
            C{i} = [ (s:e-1)' , (s+1:e)' ]; % Assign constraints
            
            % Assign node string in mesh structure
            mesh.Constraints(j).nodeStr = (s:e)';
            
            i = i + 1;
            
            j = j + 1;
            
            drawnow; % for graphics
            
        end
        
    end
    

    
    %--------------------------------------------------------
    % Internal barrier
    %--------------------------------------------------------
    
    % Find all external barrier constraints in PTS
    idx = find(ismember([PTS.Constraints.num],[4 5 24 25]))';
    
    if ~isempty(idx)
 
        for k = 1:numel(idx)
            
            % Note: weirs are not closed boundaries themselves but are
            % stored that way for computing the signed distance function.
            % For constraining we need to put the coordinates back into
            % their original format.
            
            % Store data
            mesh.Constraints(j).data = PTS.Constraints(idx(k)).data;
            
            % Store constraint value
            mesh.Constraints(j).num  = PTS.Constraints(idx(k)).num;
            
            % Store constraint value
            mesh.Constraints(j).type = PTS.Constraints(idx(k)).type;
            
            % Number of nodes in each segment
            n = (length(PTS.Constraints(idx(k)).xy(:,1)) - 1)/2;
            
            % Store first segment in pfix
            pC{i} = PTS.Constraints(idx(k)).xy(1:n,1:2);
            
            s = e + 1; % Start indexing counter
            
            % Ending point in constraints
            e = size(pC{i},1) + e;
            
            % Assign constraints
            % Do not want to constrain last point with first point
            C{i} = [ (s:e-1)' , (s+1:e)' ];
            
            % Store first node string
            nodeStr = (s:e)';
            
            % Store second segment
            i = i + 1;
                        
            % Store second segment in pfix
            pC{i} = flipud( PTS.Constraints(idx(k)).xy(n+1:end-1,1:2) );
            
            % Increase index
            s = e + 1;
            e = size(pC{i},1) + e;
            
            % Assign constraints
            % Do not want to constrain last point with first point
            C{i} = [ (s:e-1)' , (s+1:e)'];
            
            % Assign node string back to back in mesh structure
            mesh.Constraints(j).nodeStr = [nodeStr (s:e)'];
            
            i = i + 1;
            
            j = j + 1;
            
            drawnow; % for graphics

        end
                
    end
        
    % Convert cell to matrix
    C = cell2mat(C); 
        
end

% Remove duplicate constraints
C = unique(C,'rows');

%------------------------------------------------------------------------------
% Fix the end points of the open ocean boundary
%------------------------------------------------------------------------------

% Find all external barrier constraints in PTS
idx = find(ismember([PTS.Constraints.num],-1))';

if ~isempty(idx)
    
    for k = 1:noSeg
        
        % Store in pfix
        pC{i} = [PTS.Constraints(idx(k)).xy(1,:); PTS.Constraints(idx(k)).xy(end,:)];
        
        % Increase index
        s = e + 1;
        e = size(pC{i},1) + e;
        
        % Store starting and ending nodes in mesh structure.
        % The final node string will be extracted after mesh is complete.
        mesh.Constraints(j).nodeStr = (s:e)';
        
        % Store constraint value
        mesh.Constraints(j).num = -1;
        
        % Store constraint value
        mesh.Constraints(j).type = PTS.Constraints(idx(k)).type;
        
        i = i + 1;
        
        j = j + 1;
        
        drawnow; % for graphics
        
    end
    
end

%------------------------------------------------------------------------------
% Prepare outputs
%------------------------------------------------------------------------------
pC    = unique(cell2mat(pC),'rows','stable'); % Convert cell to matrix
nC    = size(pC,1);                           % Store the number of fixed nodes
if ~isempty(pC)
p     = setdiff(p,pC,'rows','stable');        % Remove duplicate nodes
p     = [pC; unique(p,'rows','stable')];      % Prepend fix points
end
end