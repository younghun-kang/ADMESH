function polygons = PTS2PointList(varargin)
% PTS2PointList - Compiles polygons in PTS data structure into a
% point in (X,Y) to the boundaries in PTS
%
% Syntax:  polygons = PTS2PointList(varargin)
%
% Inputs:
%
% Outputs:
%
% Other m-files required: none
% Subfunctions: InPolygon,bwdist,Compute_Distance_v3
% MAT-files required: none
%
% Author: Dustin West, Colton Conroy, Ethan Kubatko
% The Ohio State University
% email address: dww.425@gmail.com
% August 2013; Last revision: 08-August-2013

%------------- BEGIN CODE --------------

% Deal inputs
[PTS,delta] = deal(varargin{:});

% Index variables of contraints
if isfield(PTS,'Constraints') && ~isempty(PTS.Constraints)
    IC = find(ismember([PTS.Constraints.num],[4 5 24 25]));
    LC = find(ismember([PTS.Constraints.num],[18 19]));
else
    IC = [];
    LC = [];
end

nPTS = length(PTS.Poly); nIC = length(IC); nLC = length(LC);

nsegments = nPTS + nIC + nLC;

polygons(nsegments) = struct('xy',[],'class',[]);

%------------------------------------------------------------------------------
% Compile all polygons
%------------------------------------------------------------------------------
% PTS data structure
for i = 1:nPTS
    
    % Compute edge lengths
    eL = sqrt(diff(PTS.Poly(i).x).^2 + diff(PTS.Poly(i).y).^2);
    
    % Determine the number of points to fit between each edge
    nP = floor(eL/delta);
    
    % For any edges where nP <= 1, nP = 2. We don't need more res. here
    nP(nP<=1) = 2;
    
    % Initialize temporary cell to hold new coordinate list
    pTemp = cell(numel(nP),1);
    
    % Loope over each edge
    for j = 1:numel(nP)
        
        pTemp{j} = [...
            linspace(PTS.Poly(i).x(j), PTS.Poly(i).x(j+1),nP(j))',...
            linspace(PTS.Poly(i).y(j), PTS.Poly(i).y(j+1),nP(j))'];
        
    end
    
    % Convert cell to mat and creat unique point list
    polygons(i).xy = unique(cell2mat(pTemp),'rows','stable');
    
    % Close polygon
    polygons(i).xy(end+1,:) = polygons(i).xy(1,:);
    
    % Store class
    polygons(i).class = 'PTS';
    
end

n = nPTS;

% Barrier Constraint structure
for i = 1:nIC
    
    % Compute edge lengths
    eL = sqrt(diff(PTS.Constraints(IC(i)).xy(:,1)).^2 + ...
        diff(PTS.Constraints(IC(i)).xy(:,2)).^2);
    
    % Determine the number of points to fit between each edge
    nP = floor(eL/delta);
    
    % For any edges where nP <= 1, nP = 2. We don't need more res. here
    nP(nP<=1) = 2;
    
    n = n+1;
    
    % Initialize temporary cell to hold new coordinate list
    pTemp = cell(numel(nP),1);
    
    % Loope over each edge
    for j = 1:numel(nP)
        
        xy1 = PTS.Constraints(IC(i)).xy(j,:);
        xy2 = PTS.Constraints(IC(i)).xy(j+1,:);
        
        pTemp{j} = [...
            linspace(xy1(1), xy2(1),nP(j))',...
            linspace(xy1(2), xy2(2),nP(j))'];
        
    end
    
    % Convert cell to mat and creat unique point list
    polygons(n).xy = unique(cell2mat(pTemp),'rows','stable');
    
    % Close polygon
    polygons(n).xy(end+1,:) = polygons(n).xy(1,:);
    
    % Store class
    polygons(n).class = 'IC';
    
end

% Line Constraint structure
% Keep track of these indices. We do not need to remove interior points
for i = 1:nLC;
    
    n = n+1;
    
    polygons(n).xy = SpacePolyEqually(PTS.Constraints(LC(i)).xy,delta);
    
    polygons(n).class = 'LC';
    
end












