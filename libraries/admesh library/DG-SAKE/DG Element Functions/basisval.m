function [phi] = basisval(p,dim,X,varargin)

%--------------------------------------------------------------------------
%
%  BASIS = basisval(p,X)
%  BASIS = basisval(p,X,elem)
%
%--------------------------------------------------------------------------
%
%  Input:
%  ------
%      p: polynomial degree of the basis expansion (required).
%      X: N by d array of the N d-dimensional points at which to evaluate
%         the basis functions (required).
%   elem: element type (optional, default is the 1d line element). The ele-
%         ments in 1d, 2d, and 3d are defined as:
%           (1d) 'line', -1 <= x <= 1.
%           (2d) 'triangle', -1 <= x,y <= 1, x + y <= 0.
%           (2d) 'quadrilateral', -1 <= x,y <= 1.
%           (3d) 'tetrahedron',
%           (3d) 'pyramid',
%           (3d) 'triangular prism', -1 <= x,y,z <= 1, x + y <= 0
%           (3d) 'hexahedron' -1 <= x,y,z <= 1
%
%  Output:
%  -------
%    phi: 1 by M structure array, where M is the number of degrees of
%         freedom for an element of type elem and degree p, with fields:
%     |
%     |----- .f:    N by 1 array of the values of the basis functions
%     |             evaluated at X.
%     |
%     |----- .dfdx: N by 1 array of the values of the derivatives of the
%     |             basis functions with respect to x evaluated at X.
%     |
%     |----- .dfdy: N by 1 array of the values of the derivatives of the
%     |             basis functions with respect to y evaluated at X (only
%                   for dim > 1).
%     |
%     |----- .dfdz: N by 1 array of the values of the derivatives of the
%                   basis functions with respect to z evaluated at X (only
%                   for dim > 2).
%
%  Dependencies: principal_function => jacobi_poly
%
%--------------------------------------------------------------------------
%
%  NOTE: (1)
%        (2) The basis functions are ordered hierarchically.
%
%--------------------------------------------------------------------------
%
%  Written by Ethan Kubatko, February 8, 2012.
%
%--------------------------------------------------------------------------

% Validate data passed to function
%----------------------------------

ip = inputParser;
vp   = @(x)validateattributes(x,{'numeric'},{'scalar','integer','>=',0});
vdim = @(x)validateattributes(x,{'numeric'},{'scalar','positive','integer','<=',3});
vX   = @(x)validateattributes(x,{'numeric'},{'real','size',[NaN,dim]});
elements{1} = {'line'};
elements{2} = {'triangle','quadrilateral'};
elements{3} = {'tetrahedron','pyramid','prism','hexahedral'};
if nargin > 3
    velem = @(x)any(strcmpi(x,[elements{:}]));
    ip.addOptional('elem','line',velem);
    ip.parse(varargin{1}); elem = ip.Results.elem;
    switch elem
        case elements{1}
            x = X(:,1);
        case elements{2}
            x = X(:,1); y = X(:,2);
        case elements{3}
            x = X(:,1); y = X(:,2); z = X(:,3);
    end
else
    elem = elements{dim}(1);
end
ip.addRequired('p',vp); ip.addRequired('dim',vdim); ip.addRequired('X',vX);
ip.parse(p,dim,X);
ip.Results;

% Under construction

switch dim
    
    case 1
        %----------------------------------------------------------------------
        % 1D
        %----------------------------------------------------------------------
        
        % Initialize output
        phi(p+1) = struct('f',[],'dfdx',[]);
        
        % Compute basis
        for i = 0:p
            psi             = principal_function(i);
            phi(i+1).f      = polyval(psi,X(:,dim));
            phi(i+1).dfdx   = polyval(polyder(psi),X(:,dim));
        end
        
        
    case 2
        %----------------------------------------------------------------------
        % 2D
        %----------------------------------------------------------------------
        
        switch elem
            
            case 'triangle'
                %--------------------------------------------------------------
                % Triangle (Have to fix for singularity at y = 1);
                %--------------------------------------------------------------
                
                eta1 = (2*(1+x)./(1-y.*(y~=1))-1).*(y~=1) + (y==1)*(2*tand(45)-1);
                eta2 = y;
                
                deta1dx = 2./(1-y.*(y~=1)).*(y~=1);
                deta1dy = 2*(1+x)./(1-y.*(y~=1)).^2.*(y~=1);
                
                % Initialize index
                l = 1;
                
                for j = 0:p
                    for i = 0:j
                        psi.a = principal_function(i);     a = polyval(psi.a,eta1);
                        psi.b = principal_function(i,j-i); b = polyval(psi.b,eta2);
                        phi.f(:,l) = a.*b;
                        phi.dfdx(:,l) = polyval(polyder(psi.a),eta1).*deta1dx.*b;
                        phi.dfdy(:,l) = polyval(polyder(psi.a),eta1).*deta1dy.*b ...
                            + a.*polyval(polyder(psi.b),eta2);
                        l = l + 1;
                    end
                end
                
            case 'quadrilateral'
                %--------------------------------------------------------------
                % Quadrilateral
                %--------------------------------------------------------------
                
                % Initialize index
                l = 1;
                
                for i = 0:p
                    for j = 0:p
                        psi.a = principal_function(i); a = polyval(psi.a,x);
                        psi.b = principal_function(j); b = polyval(psi.b,y);
                        phi.f(:,l) = a.*b;
                        phi.dfdx(:,l) = polyval(polyder(psi.a),x).*b;
                        phi.dfdy(:,l) = a.*polyval(polyder(psi.b),x);
                        l = l + 1;
                    end
                end
        end
        
    case 3
        %----------------------------------------------------------------------
        % 3D.....UNDER CONSTRUCTION
        %----------------------------------------------------------------------
        %         switch elem
        %             case 'tetrahedron'
        %             case 'pyramid'
        %             case 'prism'
        %             case 'hexahedron'
        %         end
end
end



