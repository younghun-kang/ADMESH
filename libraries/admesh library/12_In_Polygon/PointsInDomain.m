function IN = PointsInDomain(X,Y,PTS)
% PointsInDomain - Determine points inside domain
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% Last revision: 25-April-2014

%---------------------------------------------------------------------
% Begin Code
%---------------------------------------------------------------------

%---------------------------------------------------------------------
% For speed up, get initial guess of points near and in domain using
% delaunay triangulation. Note: delaunayTriangulation can be used to
% determine points inside polygon but under certain circumstances the 
% isInterior does not work effectively.
%---------------------------------------------------------------------

% Total number of segments
nsegments = numel(PTS.Poly); 

% Prepare inputs for constrained delaunay triangulation
C  = cell(nsegments,1);     % Cell for constraints

% Vertically concatenate all polygons
p = [vertcat(PTS.Poly(:).x) vertcat(PTS.Poly(:).y)];

% Create unique coordinate set with back pointers.
[p,~,ic] = unique(p,'rows','stable');

% Create constraints matrix
l = 0; % Initialize cumulative sum variable
e = 0; % Initialize end index in constraints
for k = 1:nsegments
    
    l  = l + length(PTS.Poly(k).x);
    
    % Increase index
    s = e + 1; e = l;
    
    % Store constraint matrix
    C{k} = [ic(s:e-1) ic(s+1:e)];
    
end

% Convert cell to matrix. 
C = sortrows(sort(cell2mat(C),2)); 

% Check for repeating edges. If we have repeating edges then that means one
% polygon lies on another. Both repeating constraints MUST be removed from
% the constraints matrix because the "isinterior" function will not work 
% properly. 
[~,idx1]    = unique(C,'rows','first'); % Return first occurence
[~,idx2]    = unique(C,'rows','last');  % Return last occurence
C           = C(intersect(idx1,idx2),:);% Keep intersections. 

% Initialize output
IN      = false(size(X));

% Perform delaunay triangulation
dt      = delaunayTriangulation(p,C);

% Find boundaing box of domain
xmin    = min(p(:,1));
xmax    = max(p(:,1));
ymin    = min(p(:,2));
ymax    = max(p(:,2));

% Reduce computation time by finding initial points to check.
ix      = find( X >= xmin & X <= xmax & Y >= ymin & Y <= ymax );

% Find points within triangulation
ti      = pointLocation(dt,[X(ix) Y(ix)]);

% Keep only interior triangulation
inside  = find(isInterior(dt));

% Find points in domain
IN(ix)  = ismember(ti,inside);

%  figure
%  hold on
%  plot(pt(C'),pt(C'+size(pt,1)),'-r*', 'LineWidth', 2);
%  triplot(dt.ConnectivityList(inside, :),dt.Points(:,1),dt.Points(:,2));
%  plot(X(IN),Y(IN),'k.')
%  daspect([1 1 1])

end














