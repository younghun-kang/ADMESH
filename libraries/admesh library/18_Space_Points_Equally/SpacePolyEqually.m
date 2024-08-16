function p = SpacePolyEqually(points,spacing)
%------------------------------------------------------------------------------
% Space points equally
%------------------------------------------------------------------------------


% Parameterize by computing the euclidean distance between two consecutive
% points and summing the distances cumulatively.
s = [0 ; cumsum(sqrt(sum(diff(points).^2, 2)))];


% Determine the number of equally spaced points that can be fit between the
% starting and ending points with a specified element size. Need n+1 because
% we're starting from 0 when we distribute points.
MaxL = s(end);
n = ceil(MaxL/spacing);

% Distribute distance evenly from 0 to MaxL by n points.
s_spaced = linspace(0, MaxL, n)';

% Initialize new point set
% p = zeros(n, 2);
% 
% for c = 1:n
%     
%     % Determine the index of surrounding points
%     ind0 = find(s <= s_spaced(c), 1, 'last');
%     ind1 = find(s >= s_spaced(c), 1, 'first');
%     
%     if ind0 == ind1 % For beginning and ending points
%         p(c, :) = points(ind0, :);
%         continue;
%     end
%     
%     % Store the surrounding points
%     pt0 = points(ind0, :);
%     pt1 = points(ind1, :);
%     
%     % Determine the weights associated to each surrounding point
%     w0 = s_spaced(c) - s(ind0);
%     w1 = s(ind1) - s_spaced(c);
%     
%     % Compute weighted average with neighbor positions
%     p(c, :) = (pt0 * w1 + pt1 * w0) / (w0 + w1);
%     
% end

[x,y] = SpacePolyPoints(points(:,1),points(:,2),s,s_spaced); % Mex

p = [x,y];

end