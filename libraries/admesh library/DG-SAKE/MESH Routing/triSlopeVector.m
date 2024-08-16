function [ts,Sx,Sy] = triSlopeVector(MESH)

% Get the number of elements
nelems = size(MESH.ConnectivityList,1);

% Assign element coordinates
xe = reshape(MESH.Points(MESH.ConnectivityList,1),nelems,3);
ye = reshape(MESH.Points(MESH.ConnectivityList,2),nelems,3);
ze = reshape(MESH.Points(MESH.ConnectivityList,3),nelems,3);

% compute edge differences
X1 = xe(:,2) - xe(:,1); X2 = xe(:,1) - xe(:,3);
Y1 = ye(:,2) - ye(:,1); Y2 = ye(:,1) - ye(:,3);
Z1 = ze(:,2) - ze(:,1); Z2 = ze(:,1) - ze(:,3);

% Compute the steepest slope vector
i   =  (Y1.*Z2 - Y2.*Z1).*(X1.*Y2 - X2.*Y1);
j   =  (Z1.*X2 - Z2.*X1).*(X1.*Y2 - X2.*Y1);
k   = -((Z1.*X2 - Z2.*X1).^2 + (Y1.*Z2 - Y2.*Z1).^2);
ts  = [i j k];

% Normalize vector
ts = ts./repmat(sqrt(sum(ts.^2,2)),1,3);

% Compute slope in x and y direction
Sx = abs(ts(:,3)).*ts(:,1);
Sy = abs(ts(:,3)).*ts(:,2);

end
