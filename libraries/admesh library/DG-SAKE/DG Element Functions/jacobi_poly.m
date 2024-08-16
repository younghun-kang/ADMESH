function p = jacobi_poly(n,a,b)

% Syntax
%--------------------------------------------------------------------------
% p = jacobi_poly(n,a,b)
%
% Description
%--------------------------------------------------------------------------
% p = jacobi_poly(n,a,b), where n is a non-negative integer and a and b are
% real constants > -1, returns the coefficients of the n-th degree Jacobi 
% polynomial in Matlab polynomial form.
%
% Examples
%--------------------------------------------------------------------------
% Example 1: Compute the 3rd degree degree Jacobi of weights a = b = 1
%
% >> p = jacobi_poly(3,1,1)
%
% p =
%
%     7     0    -3     0
%
%--------------------------------------------------------------------------
% Example 2: Compute the 1st degree Jacobi polynomial of weights a = b = 0. 
% Note this corresponds to the 1st degree Legendre polynomial, as verified
% below
%
% >> J = jacobi_poly(1,0,0)
%
% J =
%
%     1     0
%
% >> L = legendre_poly(1)
%
% L =
%
%     1     0
%
%--------------------------------------------------------------------------
% This function was written by Ethan Kubatko
%--------------------------------------------------------------------------

% Validate input

in = inputParser;
vnab = @(x)validateattributes(x,{'numeric'},{'scalar','integer','>=',0});
in.addRequired('n',vnab); in.addRequired('a',vnab); in.addRequired('b',vnab); 
in.parse(n,a,b);
in.Results; 

% Construct the Jacobi polynomials

p = zeros([1,n+1]);
for k = 0:n
    poly = 1;
    for i = 1:k
        poly = conv(poly,[1/2,-1/2]);
    end
    p(n+1-k:n+1) = p(n+1-k:n+1) + factorial(n)/(factorial(k)*...
                        factorial(n-k))*gamma(a+b+n+k+1)/gamma(a+k+1)*poly;
end
p = gamma(a+n+1)/(gamma(a+b+n+1)*factorial(n))*p;
end