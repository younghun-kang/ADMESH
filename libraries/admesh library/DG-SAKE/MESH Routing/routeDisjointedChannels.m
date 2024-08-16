function MESH = routeDisjointedChannels(MESH)
%==========================================================================
% Route 1D cascades in mesh that do not connect to a channel or boundary.
%
% The edges identified as cascade routes will not take lateral inflow from
% neighboring elements
%
% Author:
% Dustin West
% Assistant Engineer | Hazen and Sawyer               
% dwest@hazenandsawyer.com | hazenandsawyer.com
%==========================================================================

% Assign edge type
et = MESH.edgeType;

% Assign edge connectivity
e = MESH.edges;

% Find 1D edges
loc = find(et == -2);

% If there are no 1D edges, exit the program
if isempty(loc)
    return
end

% Assign elevation
z = MESH.Points(:,3);

% Assign edges that are 1D elements
e1 = e(loc,:);

nedges = size(e1,1);

% Sort edges from low to high elevation
[~,i] = sort(z(e1),2); 

% Get the edges corresponding to the boundaries
trep    = triangulation(MESH.ConnectivityList,MESH.Points);
fb      = freeBoundary(trep);
vt      = vertexAttachments(trep); 

h = waitbar(0,'Routing 1D cascades...'); drawnow

% Loop over each edge
for k = 1:nedges
    
    if mod(k,15) == 0
        waitbar((k / nedges)); drawnow
    end
    
    % Find node with lower elevation
    nlow = e1(k,i(k,1));
    
    % Temporarily remove current edge from list
    edge = e1(k,:); e1(k,:) = [nan nan];
        
    % (1) If node is boundary node, continue to next loop.
    if any(nlow == fb(:))
        e1(k,:) = edge; % Restore edge
        continue;
    end
    
    % (2) Check if node is connected to another 1D element
    [r,c] = find(nlow == e1);
    
    if isempty(r) % No attachments, route channel
        
        e1(k,:) = edge; % Restore edge
        
        while 1
            
            % Retreive all vertice attachments
            ti = MESH.ConnectivityList(vt{nlow},:);% Elements
            va = ti(ti ~= nlow);                % vertex list
            
            % Find the node with the lowest elevation
            [~,vi] = min(z(va));
                        
            % The new 1D element edge is:
            newEdge = sort([va(vi), nlow]);
            
            % Find edge is edge list
            lie = (newEdge(1) == e(:,1) & newEdge(2) == e(:,2));
                                    
            et(lie) = -3; % Set edge flag as routing flag
            
            % Assign new node
            nlow = va(vi);
            
            % Check if new node is boundary (outflow) node
            if any(nlow == fb(:))
                break; % continue to next loop.
            end
            
            % Check if node is now connected to another 1D element
            if any(nlow == e1(:))
                break; % continue to next loop.
            end
            
        end
                
    else % Check attachment elevations in case we're at a junction
        
        % Get elevation of neighboring attachments
        c = (c + 1) - (c == 2) - (c == 2); % Logical for column switch
        zn = z(e1(r,c));
        
        % Is there any node elevation that is lower?        
        if any(zn(:) <= z(nlow))
            e1(k,:) = edge; % Restore edge
            continue
        else % Route channel
            
            e1(k,:) = edge; % Restore edge
            
            while 1
                
                % Retreive all vertice attachments
                ti = MESH.ConnectivityList(vt{nlow},:);% Elements
                va = ti(ti ~= nlow);                % vertex list
                
                % Find the node with the lowest elevation
                [~,vi] = min(z(va));
                
                % The new 1D element edge is:
                newEdge = sort([va(vi), nlow]);
                
                % Find edge is edge list
                lie = (newEdge(1) == e(:,1) & newEdge(2) == e(:,2));
                                
                et(lie) = -3; % Set edge flag as routing flag
                
                % Assign new node
                nlow = va(vi);
                
                % Check if new node is boundary (outflow) node
                if any(nlow == fb(:))
                    break; % continue to next loop.
                end
                
                % Check if node is now connected to another 1D element
                if any(nlow == e1(:))
                    break; % continue to next loop.
                end
                
            end
            
        end
        
    end
    
    % Restore edge
    e1(k,:) = edge;

end

close(h); drawnow

% Re-Save edge type ID
MESH.edgeType = et;