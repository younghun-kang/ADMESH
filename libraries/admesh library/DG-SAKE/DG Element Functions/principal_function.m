function [psi] = principal_function(p,varargin)

%--------------------------------------------------------------------------
%
%  [psi] = principal_function(p)
%  [psi] = principal_function(p,q)
%  [psi] = principal_function(p,q,r)
%   
%  Returns the (orthogonal) principal function 
%
%  psi^a_p, psi^b_pq, or psi^c_pqr,
%
%  which can be used to construct orthogonal modal basis functions for 
%  standard one-, two-, and three-dimensional element regions; see, for 
%  example, [1]. The principal functions are polynomials represented as row 
%  vectors containing the coefficients ordered by descending powers. For 
%  example, the principal function 
% 
%  psi^a_3 = 5/2*z^3 - 3/2*z 
%
%  would be represented as
%
%  psi = [5/2, 0, -3/2, 0 ].
%
%  This type of polynomial representation can be used with the various Mat-
%  lab polynomial functions, e.g., polyval, polyder, polyint, etc. 
%
%--------------------------------------------------------------------------
%
%  Input:
%  ------
%       p: degree of the one-dimensional principal function (required). 
%    q, r: indices related to the two- and three-dimensional principal
%          functions that can be used to construct orthogonal modal basis 
%          functions of degree p+q and p+q+r for standard two- and three-
%          dimensional element regions, respectively (optional).
%
%  Output:
%  -------
%     psi: 1 by N row vector of the principal function polynomial coeffici-
%          ents ordered by descending powers, where N=p+1, N=p+q+1, or 
%          N=p+q+r+1 for the 1, 2, and 3 dimensional cases, respectively.
%
%  Dependicies: jacobi_poly
%
%--------------------------------------------------------------------------
%
%  REFS: [1]  G. Karniadakis, S.J. Sherwin, Spectral/hp element methods for 
%             computational fluid dynamics, Oxford University Press, Oxford
%             (2005).
%
%--------------------------------------------------------------------------
%
%  Written by Ethan Kubatko, January 25, 2012
%
%--------------------------------------------------------------------------

% Validate data passed to function
%----------------------------------

vpqr = @(x)validateattributes(x,{'numeric'},{'scalar','integer','>=',0});
ip = inputParser;
ip.addRequired('p',vpqr); 
ip.addOptional('q',0,vpqr); ip.addOptional('r',0,vpqr);   
ip.parse(p,varargin{:});
ip.Results; q = ip.Results.q; r = ip.Results.r;
 
% Construct the principal function
%---------------------------------

psi = 1;
switch nargin
    case 1
        psi = jacobi_poly(p,0,0);
    case 2
        for n = 1:p
            psi = conv(psi,[-1/2, 1/2]);
        end
        psi = conv(psi,jacobi_poly(q,2*p+1,0));
    case 3
        for n = 1:p+q
             psi = conv(psi,[-1/2, 1/2]);
        end
        psi = conv(psi,jacobi_poly(r,2*p+2*q+2,0));
end

    