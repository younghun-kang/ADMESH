function [phi,psi,w] = DG_Element_Functions(p2,p1,MESH)
%==========================================================================
% Construct 2D and 1D basis functions.
%  NOTE: (1) The special case a = b = 0 corresponds to the n-th degree 
%            Legendre polynomial.
%
%  phi : (Basis Functions) data structure with the following fields
%        'elem2' : Basis evaluated at area integration points.
%        'edge2' : Basis evalauated at edge integration points
%        'vert2' : Basis evalauated at vertices
%
%  psi : (Transformation Functions) (1x3) data structure with the following fields
%        'elem2' : Transformation function for area integration points.
%        'edge2' : Transformation function for the edge integration points.
%        'plot' : Transformation function for the neighbor edge points
%
%    w : Weights
%
% Author:
% Dustin West
% Assistant Engineer | Hazen and Sawyer               
% dwest@hazenandsawyer.com | hazenandsawyer.com
%==========================================================================

%--------------------------------------------------------------------------
% Assemble 2D element basis and transformation functions
%--------------------------------------------------------------------------

% Get element quadrature points and corresponding weights (for triangles)
[ x.elem, y.elem, w.elem ] = triangle_quad(2*p2);

% Get edge lobatto quadrature points and corresponding weights
[xi,w.edge] = lobatto_quad(p2+2);          

% Concatenate triangular edge integration points
x.edge = [ -xi; -ones(size(xi));             xi  ];
y.edge = [  xi;             -xi; -ones(size(xi)) ];

% Create vertex points
x.vert = [-1; 1; -1]; y.vert = [-1; -1; 1];

% Construct basis functions 
phi.elem2    = basisval(p2,2,[x.elem,y.elem],'triangle');
phi.edge2    = basisval(p2,2,[x.edge,y.edge],'triangle');
phi.vert2    = basisval(p2,2,[x.vert,y.vert],'triangle');
%phi.center2  = basisval(p,2,[-1/3,-1/3],'triangle');

% Construct master element transformation functions
psi.elem2 = -1/2*[(x.elem+y.elem), -(1+x.elem), -(1+y.elem) ];
psi.edge2 = -1/2*[(x.edge+y.edge), -(1+x.edge), -(1+y.edge) ];
psi.vert2 = -1/2*[(x.vert+y.vert), -(1+x.vert), -(1+y.vert) ];

%--------------------------------------------------------------------------
% Assemble 1D element basis and transformation functions
%--------------------------------------------------------------------------
if ~isfield(MESH, 'channelNetwork'); return; end

% Get edge lobatto quadrature points and corresponding weights
%[xi,w.elem1] = lobatto_quad(p1+2); 

% Compute basis function
phi.elem1 = basisval(p1,1,xi);
phi.edge1 = basisval(p1,1,[-1;1]);

% Compute element transformation (Element points only needed)
psi.elem1 = [polyval([-1/2 1/2],xi), polyval([ 1/2 1/2],xi)];

