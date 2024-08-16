function MESH = mesh_neighbors(MESH)
%%
% Assign element edge list
E = cellfun(@(x)([x(1:end-1,1),x(2:end,1)]), MESH.channelNetwork,'uni',0);
E = vertcat(E{:});

% Get elevation
Z = reshape( MESH.Points(E',3), size(E));

% Compute the number of elements in 1D
nelems  = size(E,1); 

% Intialize neighbors matrix
N = cell(nelems,2);

% For each element...
for k = 1:nelems
   
    % Orient low to high
    if Z(k,1) > Z(k,2); E(k,:) = E(k,[2 1]); end
    
    % Assign "left" node
    node = E(k,1);
    
    % Find the neighbor(s) to the "left"
    ix      = node == E; 
    ix(k,1) = false; % Remove current find
    N{k,1}  = find( any(ix,2) );
    
    % Assign "left" node
    node = E(k,2);
    
    % Find the neighbor(s) to the "right"
    ix      = node == E;
    ix(k,2) = false; % Remove current find
    N{k,2}  = find( any(ix,2) );
    
end

MESH.Neighbors = N;

% figure; hold on
% for k = 1:nelems
%     
%     Plot element k
%     plot(X.elem1(:,k),Y.elem1(:,k),'b-.');
%     
%     Plot "left" neighbors
%     for j = 1:length(N{k,1})
%         plot(X.elem1(:,N{k,1}(j)),Y.elem1(:,N{k,1}(j)),'r-.');
%     end
% 
%     Plot "right" neighbors
%     for j = 1:length(N{k,2})
%         plot(X.elem1(:,N{k,2}(j)),Y.elem1(:,N{k,2}(j)),'g-.');
%     end
% 
%     pause;
%     cla
%     
% end


end