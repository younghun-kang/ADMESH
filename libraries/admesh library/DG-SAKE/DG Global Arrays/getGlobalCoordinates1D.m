function  [X,Y,Z,h] = getGlobalCoordinates1D(MESH,psi,ngauss,nelems)
% *************************************************************************
% Assembly/Compute the global coordinates and element sizes
%
% Author:
% Dustin West
% Assistant Engineer | Hazen and Sawyer               
% dwest@hazenandsawyer.com | hazenandsawyer.com
% *************************************************************************

% Initialize global vectors
X = zeros(nelems*(ngauss),1);
Y = zeros(nelems*(ngauss),1);
Z = zeros(nelems*(ngauss),1);

% Initialize cell for hold ing element sizes
L       = cell(size(MESH.channelNetwork));

% Initialize indexing vector
id      = 0;

% Loop over each cascade
for k = 1:length(MESH.channelNetwork)
    
    % Get (x,y) coordinates
    x = MESH.Points(MESH.channelNetwork{k}(:,1),1);
    y = MESH.Points(MESH.channelNetwork{k}(:,1),2);
    z = MESH.Points(MESH.channelNetwork{k}(:,1),3);
    
    % Get the number of points
    np = length(x)-1;
    
    % generate index vector
    id = (id(end)+1):(id(end) + ngauss*np);
    
    %     % Compute integration points
    X(id) = [psi.elem1]*[x(1:end-1,1) x(2:end,1)]';
    Y(id) = [psi.elem1]*[y(1:end-1,1) y(2:end,1)]';
    Z(id) = [psi.elem1]*[z(1:end-1,1) z(2:end,1)]';
    
    % Compute element sizes in string
    L{k} = sqrt(diff(x).^2 + diff(y).^2);
    
end

% Reshape output vectors
X = reshape(X,ngauss,nelems);
Y = reshape(Y,ngauss,nelems);
Z = reshape(Z,ngauss,nelems);

% Convert Element size cell to a vector
h = vertcat(L{:}); clear L xy



end