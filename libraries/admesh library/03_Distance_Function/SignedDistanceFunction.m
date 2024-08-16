function [D,gradD] = SignedDistanceFunction(PTS,X,Y,delta,hmax,UIFigure)
% DistanceFunction - Computes the nearest distance from a background
% point in (X,Y) to the boundaries in PTS
%
% Syntax:  [D,gradD] = SignedDistanceFunction(PTS,X,Y,delta,hmax)
%
% Inputs:
%    PTS    - ADMESH Edge Structure
%    X      - (nxm) x-coordinates to rectangular grid
%    Y      - (nxm) y-coordinates to rectangular grid
%    delta  - Grid spacing
%    hmax   - Maximum Element Size
%
% Outputs:
%    D     - (nxm) distance function (negative values indicate points
%             inside domain.
%    gradD -  Array Structure; gradient of the distance function
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

% Initialize waitbar
msg = 'Computing Distance Function.';
progdlg = uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');

% Find Points in Domain
IN = PointsInDomain3(X,Y,PTS);

% Compute nearest euclidean distance using bwdist. 
% bwdist is super fast (we're talking milliseconds. Use to compute the
% initial distance outside domain and perform computations on points within
% hmax of the domain to reduce computation time.
D = double(bwdist(IN))*delta; 

% Find all points less than hmax 
ix = find(D <= hmax);

%..........................................................................
% Younghun added
%..........................................................................
% PTS.Constraints = [];

% Space each polygon coordinate by delta and compile into data structure
p = PTS2PointList(PTS,delta);

% Find the nearest neighbor in p for each point in (X,Y)
msg = 'Computing Distance Function..';
progdlg = uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');

% Create nearest neighbor search object
ns = KDTreeSearcher(vertcat(p(:).xy),'distance','euclidean');

% Perform nearest neighbor search
[idx,D(ix)] = knnsearch(ns,[X(ix) Y(ix)]);

clear ns

nsegments = numel(p);   % Total number of segments in p
l2 = 0;                 % Initialize end index

msg = 'Computing Distance Function...';
progdlg = uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg);
progdlg.Value = 1/nsegments;

for k = 1:nsegments
    
    % Increment index vector
    l1 = l2 + 1;
    l2 = size(p(k).xy,1) + l2;
    
    % Find all common indices
    commonMembers = ismember(idx,(l1:l2)'); 
    
    % Check for matches
    if ~any(commonMembers)
        progdlg.Value = k/nsegments;
        continue; 
    end
    
    % Store indexes to p(k).xy, Rescale 1:numel(p(k).xy)
    id = idx(commonMembers); id = id - (min(id)-1);
    
    % Stor indices for [X,Y]
    nn = ix(commonMembers);

    % Compute actual distance to edge segments
    D(nn) = min(distFunSubroutine(X(nn),Y(nn),p(k).xy(:,1),p(k).xy(:,2),id),D(nn));
    
    progdlg.Value = k/nsegments;
    
end
close(progdlg);
drawnow;
% Create signed distance function
D(IN) = -D(IN);

%..Younghun: Overwrite signed distance function
if length(PTS.Poly) == 1
    XYb1 = {[PTS.Poly.x(:), PTS.Poly.y(:)]};
else
    XYb1 = struct2table(PTS.Poly,'AsArray',1);
    XYb1 = cellfun(@(x,y) [x, y],XYb1.x,XYb1.y,'UniformOutput',0);
end
if isempty(PTS.Constraints)
    XYb2 = [];
else
    XYb2 = struct2table(PTS.Constraints,'AsArray',1);
    XYb2 = XYb2.xy;
end
% [~,~,D] = VectorDistanceTransform([XYq1;XYq2],[X(:),Y(:)]);
[~,~,D] = Compute8SSED_v3([XYb1;XYb2],X(:),Y(:),delta);

D = reshape(D,size(X));
D(IN) = -D(IN);



%------------------------------------------------------------------------------
% Compute the gradient
%------------------------------------------------------------------------------
[LY,LX] = size(D);

gradD.x = zeros(LY,LX);
gradD.y = zeros(LY,LX);

gradD.x(3:LY-2,3:LX-2) = ( 1*D(3:LY-2,1:LX-4) + -8*D(3:LY-2,2:LX-3) + ...
    8*D(3:LY-2,4:LX-1) + -1*D(3:LY-2,5:LX) )/(12*delta);
gradD.y(3:LY-2,3:LX-2) = ( 1*D(1:LY-4,3:LX-2) + -8*D(2:LY-3,3:LX-2) + ...
    8*D(4:LY-1,3:LX-2) + -1*D(5:LY,3:LX-2) )/(12*delta);
