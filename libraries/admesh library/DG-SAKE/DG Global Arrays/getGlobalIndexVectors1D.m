function [IL,IR,BC] = getGlobalIndexVectors1D(MESH,npoints,ngauss)
%
% Author:
% Dustin West
% Assistant Engineer | Hazen and Sawyer               
% dwest@hazenandsawyer.com | hazenandsawyer.com

%------------------------------------------------------------------
% Initialize variables
%------------------------------------------------------------------
IL    = cell(npoints,1);  % Cell for indexing of left elements
IR    = cell(npoints,1);  % Cell for indexing of right elements
BC    = -ones(npoints,1); % Vector for boundary conditions
lv    = 1;                % Global 1D index for "left"
rv    = 1;                % Global 1D index for "right"

%------------------------------------------------------------------
% Initialize left and right indexing ( not including junctions )
%------------------------------------------------------------------
for k = 1:length(MESH.channelNetwork)
    
    % Get the number of points in string
    n = size(MESH.channelNetwork{k},1);
    
    % Populate left and right global index vectors
    for j = 1:n
        
        % Get index into IL & IR & boundary condition
        i     = MESH.channelNetwork{k}(j,2);
        BC(i) = MESH.channelNetwork{k}(j,4);
        
        
        % Assign the "left" and "right" indices based on the node we are at
        
        if (BC(i) == 2) && (j == n)
            
            % At a juncton (BC == 2) at the end of the node string (j == n)
            % This is an outflow node. Store "left" index
            
            % Assign left index value, append to existing values
            IL{i} = [IL{i} lv];
                        
            % Index lv only by 1. The "right" element is assigned later in
            % the loop under the condition (BC == 2) && (j == 1).
            lv = lv + 1;
            continue;
            
        elseif (BC(i) == 1) && (j == n)
            
            % At an outflow node (BC == 1) at the end of the node
            % string (j == n). This is an outflow node.
            % Store "left" index and copy to "right" index.
            
            % Assign index value, append to existing values
            IL{i} = [IL{i} lv];
            IR{i} = [IR{i} lv];
            
            % Incremement left global index value
            lv = lv + 1;
            continue;
            
        elseif (BC(i) == 2) && j == 1
            
            % If the first nodes' boundary condition is 2 then we are at a
            % junction and this node is going to be the node that takes in
            % the flow from the other nodes that drain into it. This index
            % is going to be stored as a "right" element index.
            
            % Assign right index value only, append to existing values
            IR{i} = [IR{i} rv];
            
        else % This is the case of BC == -1 || BC == 0
            
            % Assign index value, append to existing values
            IL{i} = [IL{i} lv];
            IR{i} = [IR{i} rv];
            
        end
        
        % Incremement left global index value
        lv = lv + ngauss - (j == 1);
        rv = rv + ngauss;
        
        
    end

end

%------------------------------------------------------------------
% Now loop over each reach in MESH.cascadeNetwork again and included
% junctions
%------------------------------------------------------------------
for k = 1:length(MESH.channelNetwork)
    
    % Get node string bc
    bc = MESH.channelNetwork{k}(:,4);
    
    % Is there a junction at the outflow?
    if bc(end) == 2
        
        % Assign local node number
        ln = MESH.channelNetwork{k}(end,2);
        
        % Get the junction index
        junc = MESH.channelNetwork{k}(end,1);
        
        % Loop over each cascade and find the corresponding inflow node
        for i = 1:length(MESH.channelNetwork)
            
            if i == k; continue; end % Skip current node string
            
            % When we have a match
            if MESH.channelNetwork{i}(1,1) == junc
                
                % Get the "right" local node number
                rn = MESH.channelNetwork{i}(1,2);
                
                % Assign value into "right" global index
                IR{ln} = [IR{ln} IR{rn}];
                
                % Since we're only looking for one node we can exit the for
                % loop
                break;
                
            end
            
        end
    end
    
    % Is there a junction at the inlet?
    if bc(1) == 2
        
        % Assign local node number
        ln = MESH.channelNetwork{k}(1,2);
        
        % Get the junction index
        junc = MESH.channelNetwork{k}(1,1);
        
        % Loop over each cascade and find the corresponding outflow nodes
        for i = 1:length(MESH.channelNetwork)
            
            if i == k; continue; end % Skip current node string
            
            % When we have a match
            if MESH.channelNetwork{i}(end,1) == junc
                
                % Get the "left" local node number
                rn = MESH.channelNetwork{i}(end,2);
                
                % Assign value into "right" global index
                IL{ln} = [IL{ln} IL{rn}];
                
            end
            
        end
        
    end

end

end