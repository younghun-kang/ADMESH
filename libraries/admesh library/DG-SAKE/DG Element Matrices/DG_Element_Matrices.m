function [a,b,c] = DG_Element_Matrices(phi,w)
%==========================================================================
% Construct DG element matrices
%
% Author:
% Dustin West
% Assistant Engineer | Hazen and Sawyer               
% dwest@hazenandsawyer.com | hazenandsawyer.com
%==========================================================================

%--------------------------------------------------------------------------
% Construct 2D element matrices
%--------------------------------------------------------------------------

% Compute 2D mass matrix 
m = ( phi.elem2.f' )*( diag(w.elem) )*( phi.elem2.f );

% Compute element advection matrix
a.x = m\( [phi.elem2.dfdx]'*diag(w.elem) );
a.y = m\( [phi.elem2.dfdy]'*diag(w.elem) );

% Compute element boundary matrices
n = length(w.edge); i = 1:n:3*n;
b.ii{1} = -m\( phi.edge2.f(i(1):i(1)+n-1,:)'*diag(w.edge) ); % edge 1
b.ii{2} = -m\( phi.edge2.f(i(2):i(2)+n-1,:)'*diag(w.edge) ); % edge 2
b.ii{3} = -m\( phi.edge2.f(i(3):i(3)+n-1,:)'*diag(w.edge) ); % edge 3

% Compute element source matrix
c.ii = m\( phi.elem2.f'*diag(w.elem) );

%--------------------------------------------------------------------------
% Construct 1D element matrices
%--------------------------------------------------------------------------
if ~isfield(phi, 'elem1'); return; end

% Compute 1D mass matrix
m = ( [phi.elem1.f]' )*( diag(w.edge) )*( [phi.elem1.f] );

% Compute the advection matrix
a.i = m\( [phi.elem1.dfdx]'*diag(w.edge) );

% Compute boundary vectors
b.i = m\( [phi.edge1.f]' ); b.i(:,2) = -b.i(:,2);

% Compute source matrix
c.i = m\( [phi.elem1.f]'*diag(w.edge) );