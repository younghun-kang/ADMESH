function [A,B,C,IL,IR,BC,X,Y,Z,PHI,N,ZL,ZR] = DG_Global_Assembly(a,b,c,phi,psi,MESH)
%==========================================================================
% Construct DG global arrays
%
% Author:
% Dustin West
% Assistant Engineer | Hazen and Sawyer               
% dwest@hazenandsawyer.com | hazenandsawyer.com
%==========================================================================

%==========================================================================
% Construct 2D global arrays
%==========================================================================

% Initialize variables
nbp         = size(b.ii{1},2);                  % # of boundary points
ndof        = size(phi.elem2.f,2);              % # of degrees of freedom
nelems      = size(MESH.ConnectivityList,1);    % # of elements

% Assign element coordinates
xe = reshape(MESH.Points(MESH.ConnectivityList,1),nelems,3)';
ye = reshape(MESH.Points(MESH.ConnectivityList,2),nelems,3)';
ze = reshape(MESH.Points(MESH.ConnectivityList,3),nelems,3)';

% Compute edge lengths
L{1} = sqrt((xe(3,:)-xe(2,:)).^2 + (ye(3,:)-ye(2,:)).^2);
L{2} = sqrt((xe(1,:)-xe(3,:)).^2 + (ye(1,:)-ye(3,:)).^2);
L{3} = sqrt((xe(2,:)-xe(1,:)).^2 + (ye(2,:)-ye(1,:)).^2);

% Compute element area
Ae      = polyarea(xe,ye);

%--------------------------------------------------------------------------
% Compute element Jacobians
%--------------------------------------------------------------------------
J(1).x = spdiags(((ye(3,:)-ye(1,:) )./(Ae))',0,nelems,nelems);
J(2).x = spdiags(((ye(1,:)-ye(2,:) )./(Ae))',0,nelems,nelems);
J(1).y = spdiags(((xe(1,:)-xe(3,:) )./(Ae))',0,nelems,nelems);
J(2).y = spdiags(((xe(2,:)-xe(1,:) )./(Ae))',0,nelems,nelems);

%--------------------------------------------------------------------------
% Compute Global Advection Matrices
%--------------------------------------------------------------------------
A.x = kron(J(1).x,sparse(a.x)) + kron(J(2).x,sparse(a.y));
A.y = kron(J(1).y,sparse(a.x)) + kron(J(2).y,sparse(a.y));

clear J

%--------------------------------------------------------------------------
% Compute Global Source Matrix
%--------------------------------------------------------------------------
C.ii = kron(spdiags(ones(nelems,1),0,nelems,nelems),c.ii);

%--------------------------------------------------------------------------
% Compute Global Boundary Arrays
%--------------------------------------------------------------------------
[B,IL,IR,BC] = assembleBoundaryMat2D(MESH,b,Ae,L,ndof,nbp,nelems);

%--------------------------------------------------------------------------
% Compute Global Normal Arrays
%--------------------------------------------------------------------------

% Compute edge normals
nx1 = (ye(3,:)-ye(2,:))./L{1}; ny1 = (xe(2,:)-xe(3,:))./L{1};
nx2 = (ye(1,:)-ye(3,:))./L{2}; ny2 = (xe(3,:)-xe(1,:))./L{2};
nx3 = (ye(2,:)-ye(1,:))./L{3}; ny3 = (xe(1,:)-xe(2,:))./L{3};

% Replicate values for each quadrature point per edge
N.x = reshape([ones(nbp,1)*nx1;ones(nbp,1)*nx2;ones(nbp,1)*nx3],nelems*nbp*3,1);
N.y = reshape([ones(nbp,1)*ny1;ones(nbp,1)*ny2;ones(nbp,1)*ny3],nelems*nbp*3,1);

clear nx1 nx2 nx3 ny1 ny2 ny3

%--------------------------------------------------------------------------
% Assemble 2D Global Element Coordinates
%--------------------------------------------------------------------------
X.elem2 = psi.elem2*xe; Y.elem2 = psi.elem2*ye; Z.elem2 = psi.elem2*ze;
X.edge2 = psi.edge2*xe; Y.edge2 = psi.edge2*ye; Z.edge2 = psi.edge2*ze;
X.vert2 = psi.vert2*xe; Y.vert2 = psi.vert2*ye; Z.vert2 = psi.vert2*ze;

%--------------------------------------------------------------------------
% Assemble the block diagonal global PHI matrices
%--------------------------------------------------------------------------
PHI.elem2 = kron(spdiags(ones(nelems,1),0,nelems,nelems),phi.elem2.f);
PHI.edge2 = kron(spdiags(ones(nelems,1),0,nelems,nelems),phi.edge2.f);
PHI.vert2 = kron(spdiags(ones(nelems,1),0,nelems,nelems),phi.vert2.f);
%PHI.center2  = kron(spdiags(ones(nelems,1),0,nelems,nelems),phi.center2.f);

%------------------------------------------------------------------
% Assemble the block diagonal global PSI matrices
%------------------------------------------------------------------
PSI.elem2 = kron(spdiags(ones(nelems,1),0,nelems,nelems),psi.elem2);
PSI.edge2 = kron(spdiags(ones(nelems,1),0,nelems,nelems),psi.edge2);

%==========================================================================
% Construct 1D global arrays
%==========================================================================
if ~isfield(MESH,'channelNetwork'); [ZL,ZR] = deal(double.empty); return; end

% Initialize variables
npoints = sum(cellfun(@(x) size(x,1),MESH.channelNetwork));     % # of points
nelems  = sum(cellfun(@(x) size(x,1),MESH.channelNetwork)-1);   % # of elements
ngauss  = size(a.i,2);                                          % # of gauss points
ndof    = size(a.i,1);                                          % # of degrees of freedom

%--------------------------------------------------------------------------
% Get global 1D coordinates and element sizes
%--------------------------------------------------------------------------
[X.elem1,Y.elem1,Z.elem1,h] = getGlobalCoordinates1D(MESH,psi,ngauss,nelems);

%--------------------------------------------------------------------------
% Compute element Jacobians
%--------------------------------------------------------------------------
J = spdiags((2./h),0,nelems,nelems);

%--------------------------------------------------------------------------
% Compute Global Advection Matrix
%--------------------------------------------------------------------------
A.i = kron( J , sparse(a.i) );

%--------------------------------------------------------------------------
% Compute Global Source Matrix
%--------------------------------------------------------------------------
C.i = kron( spdiags(ones(nelems,1),0,nelems,nelems) , sparse(c.i) );

%------------------------------------------------------------------
% Assemble left & right gauss end points & BC
%------------------------------------------------------------------
[IL.i,IR.i,BC.i] = getGlobalIndexVectors1D(MESH,npoints,ngauss);

%------------------------------------------------------------------
% Assemble Boundary Flux Matrix
%------------------------------------------------------------------
[B.i,IL,IR,ZL,ZR] = assembleBoundaryMat1D(MESH,IL,IR,Z,b,h,ndof,npoints,ngauss,nelems);

%------------------------------------------------------------------
% Create the global PHI matrix for the elements
%------------------------------------------------------------------
PHI.elem1 = kron(spdiags(ones(nelems,1),0,nelems,nelems),[phi.elem1.f]);

%------------------------------------------------------------------
% Create the global PHI matrix for the elements
%------------------------------------------------------------------
PHI.edge1 = kron(spdiags(ones(nelems,1),0,nelems,nelems),[phi.edge1.f]);
