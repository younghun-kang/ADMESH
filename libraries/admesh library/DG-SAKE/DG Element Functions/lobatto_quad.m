function [x,w] = lobatto_quad(n)

%--------------------------------------------------------------------------
%
%  [x,w] = lobatto_quad(n)
%   
%  Returns the n Lobatto quadrature points, x (in ascending order) and
%  their respective weights w to approximate an integral of the form
%
%  int( f(x),x,-1,1 ). 
%
%--------------------------------------------------------------------------
%
%  Input:
%  ------
%    n: number of points.
%
%  Output:
%  -------
%    x: 1 by n row vector of the n Lobatto integration points
%    w: 1 by n row vector of the respective weights
%
%
%  Dependencies: jacobi_poly
%
%--------------------------------------------------------------------------
%  NOTES:  (1) An n-point Lobatto quadrature will integrate all
%              polynomials f(x) up to (and including) degree 2*n-3 exactly.
%          (2) The integration points include the endpoints -1 and 1.
%--------------------------------------------------------------------------
%
%  Written by Ethan Kubatko
%
%--------------------------------------------------------------------------

% Validate data passed to function
%----------------------------------

vn = @(x)validateattributes(x,{'numeric'},{'scalar','integer','>',1});
ip = inputParser;
ip.addRequired('n',vn); 
ip.parse(n);
ip.Results;

% Construct the (n-1)-th degree Jacobi polynomial derivative
%-----------------------------------------------------------

P = jacobi_poly(n-1,0,0); dPdx = polyder(P);

% Obtain the n Lobatto integration points and sort in ascending order
%--------------------------------------------------------------------

x = sort([-1;roots(dPdx);1]);

% Compute the corresponding n weights
%------------------------------------

w = [2/(n*(n-1)); 2./(n*(n-1)*polyval(P,x(2:n-1)).^2); 2/(n*(n-1))];
end