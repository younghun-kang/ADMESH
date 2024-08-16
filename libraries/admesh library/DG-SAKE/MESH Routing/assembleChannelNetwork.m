function MESH = assembleChannelNetwork(MESH)
%==========================================================================
% Assemble cascade network from the edges in MESH identified as cascade
% edges. 
%
% The output will be 2 additional fields add to MESH
%
%   (a) cascadeNetwork: Nx1 cell containing Mx2 vectors.
%                       Where N is the total numer if cascade strings.
%                       M is the number of points in each cascade string.
%                       The 1st column in each string references the index
%                       in MESH.Points for the coordinates. The 2nd column 
%                       is used for the global assembly of the cascade 
%                       network. The 3rd column will have values ranging
%                       from 1 to NP, where NP is the total number of
%                       points in the channel network.
%
%   (b) edge2cascade: Nx1 vector (where N is to total number of edges in
%                     MESH) containing reference the the cell in
%                     cascadeNetwork with the scascade string that contains
%                     that edge. 
%                       
%
% Author:
% Dustin West
% Assistant Engineer | Hazen and Sawyer               
% dwest@hazenandsawyer.com | hazenandsawyer.com
%==========================================================================

%==========================================================================
% Check for 1D channels
%==========================================================================
et      = MESH.edgeType;                % Assign edge type
loc     = find(et == -2 | et == -3);    % Find 1D edges

% If there are no 1D edges, exit the program
if isempty(loc)
    return
end

%==========================================================================
% Initialize variables
%==========================================================================

z       = MESH.Points(:,3);             % Assign elevation
e       = MESH.edges(loc,:);            % Assign 1D edges, sort columns
[ix,ij] = getBranchPoints(e);           % Find end points and junction points
cNetwork = cell(numel(ix),1);           % Initialize cell for holding channels
edge2cascade = zeros(size(et));         % Global edge 2 cascade reference

%==========================================================================
% Group channels based on edge index
%==========================================================================
j = 1; % Index for channel network
repeatFlag = nan; count = 0;

fprintf('Assembling 1D channel network...\n\n');
h = waitbar(0,'Assembling 1D channel network...'); drawnow

% We have an outer while loop because we'll want to loop through ix several
% times
while ~all(isnan(e(:)))
    
    if sum(isnan(e(:))) == repeatFlag
        repeatFlag = sum(isnan(e(:))); 
        count = count + 1;
        if count >= 10; break; end
    else
        repeatFlag = sum(isnan(e(:))) ; count = 0;
    end
        
    waitbar((1 / length(ix)),h,...
        ['Assembling 1D channel network... ' num2str(sum(isnan(e(:)))./numel(e),'%.2f'), '%'])
    
    for k = 1:length(ix)
        
        % Update waitbar
        if mod(k,15) == 0
            waitbar((k / length(ix))); drawnow
        end

        % Initialize edge index (ei) as empty
        str = zeros(length(loc),1); n = 1;
        
        % Find the edge associated with node ix(k)
        [r,c] = find(ix(k) == e);

        % Continue on, node already taken care of
        if isempty(r); continue; end
        
        % Assign node
        str(n) = e(r(1),c(1)); n = n + 1;
        
        % Assemble channel
        while 1
            
            %------------------------------------------------------------
            % Mark edge location in channel for lateral inflow
            %------------------------------------------------------------
            if et(loc(r(1))) == -2
                edge2cascade(loc(r(1))) = j;
            end
                        
            % Update column index
            c = (c(1) + 1) - (c(1) == 2) - (c(1) == 2); % Logical for column switch
            
            % Get next node to look for
            nn = e(r(1),c(1));
            
            %------------------------------------------------------------
            % Perform a check in the case where we reach a ridge
            %------------------------------------------------------------
            if n >= 3
                
                % Compute the sign of the slope of the last segment and
                % compare the sign of this slope with the sign of the slope
                % about the be appended
                
                % If it is zero, cut of the node string here
                if sum(sign(diff(z([str((n-2):(n-1));nn])))) == 0
                    
                    % Remove trailing zeros
                    str = nonzeros(str);
                    
                    % Flip string?
                    % NOTE should be (<) but the = was added for v-shaped test case
                    if z(str(1)) < z(str(end))
                        %disp(abs(z(str(1)) - z(str(end))))
                        str = flipud(str);
                    end
                    
                    % Store channel network
                    cNetwork{j} = nonzeros(str); j = j + 1;
                    
                    % Continue to next loop
                    break
                    
                else
                    % Store node
                    str(n) = nn; n = n + 1;
                end
                
            else
                
                % Store node
                str(n) = nn; n = n + 1;
                
            end
                   
            %------------------------------------------------------------
            % Remove row from e matrix
            %------------------------------------------------------------
            e(r(1),:) = nan;
            
            % Find the edge associated with node ix(k)
            [r,c] = find(str(n-1) == e);
            
            % Conditional for moving to next channel
            if isempty(r) || numel(r) > 1 || any(e(r(1),c(1)) == ij)
                
                % Remove trailing zeros
                str = nonzeros(str);
                
                % Flip string?
                % NOTE should be (<) but the = was added for v-shaped test case 
                if z(str(1)) <= z(str(end))
                    %disp(abs(z(str(1)) - z(str(end))))
                    str = flipud(str);
                end
                
                % Store channel network
                cNetwork{j} = nonzeros(str); j = j + 1;
                
                break
            end
            
        end
        
    end

end

close(h); drawnow

% Remove empty cells
i = cellfun(@isempty,cNetwork);
cNetwork(i) = [];

% Set routing edges that take no lateral inflow (et == -3) back to et ==
% -1. This edge is still a flow edge but for the purpose of routing 1D
% channels was flagged as -3.
%et(et == -3) = -1; 

%==========================================================================
% Append nodal numbering system & BC for global assembly
%==========================================================================

% Create a single array compiling all strings in cNetwork
globalIndex = vertcat(cNetwork{:});

% Use the unique function to obtain a unique indexing array starting from 1
[~,~,IC] = unique(globalIndex,'stable');

% Append vector to each cell
j = 1;
for k = 1:length(cNetwork)
    
    % Get the length of the current vector
    L = length(cNetwork{k}) + j - 1;
    
    % Determine boundary conditions
    %   -1: Flow
    %    0: No Flow
    %    1: Radiation
    %    2: Junction
    BC      = -ones(length(j:L),1);
    BC(1)   = any(cNetwork{k}(1) == ij)*2;    % Can be no flow or junction
    BC(end) = any(cNetwork{k}(end) == ij)*1+1;% Can be radiation or junction

    % Append
    %cNetwork{k} = [cNetwork{k} IC(j:L) BC];
    % Column 1: Index values to MESH.Points
    % Column 2: Global 1D node number
    % Column 3: Global 1D node numbers with no repititions. 
    % Column 4: Boundary condition
    cNetwork{k} = [cNetwork{k} (j:L)' IC(j:L) BC];
    
    % increment j
    j = j + size(cNetwork{k},1);

end

%==========================================================================
% Append channel properties to MESH
%==========================================================================
MESH.edgeType       = et;
MESH.channelNetwork = cNetwork;
%MESH.edge2cascade   = edge2cascade;

end


%==========================================================================
% SUBROUTINE
%==========================================================================
function [ix,ij] = getBranchPoints(E)
%==========================================================================
% Get branch points
%
% Author: Dustin West
%==========================================================================

E = E(:);

% Get vector of non-repeating elements and index vector
C = unique(E,'first');

count = zeros(size(C));


for k = 1:length(C)
    count(k) = sum(C(k) == E);
end


% Find indices that do not repeat (end points of segments)
ix = C(count == 1);

% Find indices where two streams exit the domain at the same point

% Find indices that are junctions
ij = C(count >= 3);

% Append ij to the bottom of ix.
ix = [ix;ij];

end