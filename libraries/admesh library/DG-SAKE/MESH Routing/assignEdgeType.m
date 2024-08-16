function MESH = assignEdgeType(MESH,tol)
%==========================================================================
% Classify edges in MESH based on their flow type.
%
%   Edge types: (a) Flow edge               ; type number: -1
%               (b) No flow (ridge) edge    ; type number:  0
%               (c) Radiation edge          ; type number:  1
%               (d) channel edge            ; type number: -2
%
% The edge types identified will be routed using the routeCascade algorithm
% and used for transporting flow in the domain where there is converging
% flow in neighboring elements.
%
% Input:    MESH - structure array with fields:
%               Points
%               ConnectivityList
%               Constraints         - struct array with fields:
%                   xy      - coordinates of constraint
%                   num     - number designation
%                   type    - string description of constraint
%                   data    - additional data pertaining to each node
%
% Output:   MESH - (1xn) structure array with fields: (n denotes number of MESHES)
%               Points
%               ConnectivityList
%               Dimension           - 1, 2, or 3
%               Constraints         - struct array with fields:
%                   xy      - coordinates of constraint
%                   num     - number designation
%                   type    - string description of constraint
%                   data    - additional data pertaining to each node
%               edges - Edge List
%               edgeAttachments - Edge Attachments (ordered outflow -> inflow)
%               edgeType - Boundary condition classification
%                            0 : No Flow (ridges)
%                            1 : Radiation edge
%                           -1 : Flow Boundary
%                           -2 : Cascade edge
%
% Author:
% Dustin West
% Assistant Engineer | Hazen and Sawyer               
% dwest@hazenandsawyer.com | hazenandsawyer.com
%==========================================================================

fprintf('Classifying element edges...\n\n');

%--------------------------------------------------------------------------
% Get geometric properties of mesh
%--------------------------------------------------------------------------
trep        = triangulation(MESH.ConnectivityList,MESH.Points);
E           = edges(trep);                             % Edge list
EA          = edgeAttachments(trep,E);                 % Edge attachements
[~,Sx,Sy]   = triSlopeVector(MESH);                    % Element slope vector

%--------------------------------------------------------------------------
% Construct edge type vector (Order EA outflow -> inflow)
%--------------------------------------------------------------------------
ET          = -ones(size(E,1),1);   % Initialize edge type vector
flowMap     = cell(size(trep,1),1);
 
%--------------------------------------------------------------------------
% Flag channels provided in MESH.Constraints
%--------------------------------------------------------------------------
% Get channel edge list
CE = cell(length(MESH.Constraints),1);
for k = 1:length(CE)
    if MESH.Constraints(k).num == 18
       CE{k} =  [...
           MESH.Constraints(k).nodes(1:end-1) ...
           MESH.Constraints(k).nodes(2:end)];
    end 
end
CE      = vertcat(CE{:}); % Channel edge list
if ~isempty(CE)
idx     = ismember(sort(E,2),sort(CE,2),'rows'); % Find edges in list
ET(idx) = -2; % Flag
end
clear idx CE

%--------------------------------------------------------------------------
% Assign coordinates
%--------------------------------------------------------------------------

% Assign x edge nodes ( for parallel computing )
x1 = MESH.Points(E(:,1),1);
x2 = MESH.Points(E(:,2),1);

% Assign y-nodes ( for parallel computing )
y1 = MESH.Points(E(:,1),2);
y2 = MESH.Points(E(:,2),2);

% Assign the 3 vertices for each element
p1 = MESH.Points(MESH.ConnectivityList(:,1),[1 2]);
p2 = MESH.Points(MESH.ConnectivityList(:,2),[1 2]);
p3 = MESH.Points(MESH.ConnectivityList(:,3),[1 2]);

for k = 1:size(EA,1)     % Loop over each row in EA
    
    switch length(EA{k})
        
        case 1 % Free boundary edge
            
            % Retrieve element indexing
            iL = EA{k};
            
            % Compute [mid-point + slope component] of the edge
            mx = (x1(k) + x2(k))/2 + Sx(iL);  %#ok<*PFBNS>
            my = (y1(k) + y2(k))/2 + Sy(iL);
            
            % Test if [mid-point + slope component] is in "left" triangle
            [IN,ON] = pointInTri([mx, my],p1(iL,:),p2(iL,:),p3(iL,:));
            
            % For boundaries, if vector is pointing into the triangle, it
            % is considering no-flow. If vector is pointing out of the
            % triangle, it is considered radiation.
            ET(k) = ~(IN || ON);
            
        case 2 % Interior edge
            
            % Retrieve "left" & "right" element indexing
            iL = EA{k}(1); iR = EA{k}(2);
            
            % Compute mid-point of the edge
            mx = (x1(k) + x2(k))/2; my = (y1(k) + y2(k))/2;
            
            % Test if [mid-point + slope component] is in "left" triangle
            [inL,onL] = pointInTri([mx + Sx(iL), my + Sy(iL)], p1(iL,:),p2(iL,:),p3(iL,:));
            
            % Test if [mid-point + slope component] is in "right" triangle
            [inR,onR] = pointInTri([mx + Sx(iR), my + Sy(iR)], p1(iR,:),p2(iR,:),p3(iR,:));
            
            if ~inL && ~onL && ~inR && ~onR     % Converging edge (1D element)
                
                % Check difference between points
                val = max(abs(...
                    abs([mx + Sx(iL), my + Sy(iL)]) - ...
                    abs([mx + Sx(iR), my + Sy(iR)])));
                if val <= tol
                    ET(k) = -1; % Edge type flag
                else
                    ET(k) = -2; % Assign channel flag
                    flowMap{iL} = [flowMap{iL} iR]; % Add contributing elem
                    flowMap{iR} = [flowMap{iR} iL]; % Add contributing elem
                end
                
                continue;
                
            elseif ~inL && onL && ~inR && ~onR  % Flow to no flow: case 1

                ET(k) = -2; % Edge type flag
                EA{k} = [iR iL]; % Order from outflow -> inflow
                flowMap{iL} = [flowMap{iL} iR]; % Add contributing elem
                continue;
                
            elseif ~inL && ~onL && ~inR && onR  % Flow to no flow: case 2

                ET(k) = -2; % Edge type flag
                flowMap{iR} = [flowMap{iR} iL]; % Add contributing elem
                continue;
                
            elseif inL && ~onL && inR && ~onR   % Ridge edge
                
                ET(k) = 0; % Edge type flag
                continue;
                
            elseif inL && ~inR                  % Flow edge

                EA{k} = [iR iL]; % Order from outflow -> inflow
                flowMap{iL} = [flowMap{iL} iR]; % Add contributing elem
            else
                flowMap{iR} = [flowMap{iR} iL]; % Add contributing elem
            end
            
    end
    
end

% Save edges
MESH.edges = E;

% Save edge attachments
MESH.edgeAttachments = EA;

% Save edge type ID
MESH.edgeType = ET;

% Save flow map for MESH
MESH.flowMap = flowMap;


end
