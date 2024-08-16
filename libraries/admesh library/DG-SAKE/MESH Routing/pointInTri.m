function [IN,ON] = pointInTri(pt,v1,v2,v3)
%==========================================================================
% Determines whether point pt is inside a triangle with vertices v1, v2, v3
%
% Author: Dustin West
%==========================================================================

% Move origin to one of the vertices
v2 = v2 - v1;
v3 = v3 - v1;
pt = pt - v1;

% Compute denominator for calculating barycentric weights
d = v2(1)*v3(2) - v3(1)*v2(2);

% Calculating barycentric weights
w1 = (pt(1)*(v2(2)-v3(2)) + pt(2)*(v3(1)-v2(1)) + v2(1)*v3(2)-v3(1)*v2(2))/d;
w2 = (pt(1)*v3(2) - pt(2)*v3(1))/d;
w3 = (pt(2)*v2(1) - pt(1)*v2(2))/d;

% Cap to zero
if abs(w1) < eps; w1 = 0; end 
if abs(w2) < eps; w2 = 0; end 
if abs(w3) < eps; w3 = 0; end 

% The point is IN if weights are (0,1)
IN = all([w1 w2 w3] > 0 & [w1 w2 w3] < 1);

% The point is ON if weights are (0,1) and one point is 0
ON = any([w1 w2 w3] == 0) && all([w1 w2 w3] >= 0) && all([w1 w2 w3] <= 1);
