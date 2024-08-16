function [fzL,fzR] = InverseCSTopo(fz,MaxChannelWidth,CoordinateSystem)
% Compute inverse of cross-section topography
% Outputs are functions with input of height and output of cross sectional 
% coordinates for left and right side, respectively

% Set arbitrary cross-sectional interval
d = linspace(-MaxChannelWidth,MaxChannelWidth,1e3);
if strcmpi(CoordinateSystem,'Unprojected (decimal degree)')
    d = d*km2deg(1e-3);
end
d = d(:);

% Set center of cross-section
I = find(abs(d) < 10);
[~, id_center] = min(fz(d(I)));
I = I(id_center);
csd_center = d(I);

% Find left/right heights with filling sinks so that inverse function can
% be defined
d_left = d(d <= csd_center);
d_left = flip(d_left);
H_left = zeros(size(d_left));
for i = 1 : length(d_left)
    H_left(i) = max(fz(d_left(1:i)));
end

d_right = d(d >= csd_center);
H_right = zeros(size(d_right));
for i = 1 : length(d_right)
    H_right(i) = max(fz(d_right(1:i)));
end

% Remove points with same height (remaining one point closest to center)
I = find(diff(H_left) > 0)+1;
fzL = @(h) interp1(H_left(I),d_left(I),h);

I = find(diff(H_right) > 0)+1;
fzR = @(h) interp1(H_right(I),d_right(I),h);

