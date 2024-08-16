function plist = Mesh2PolyList(tri,p)
% Mesh2PolyList - Segments boundaries in mesh into a NaN seperated 
% point list
%
% Syntax:  plist = Mesh2PolyList(tri,p)
%
% Inputs:
%    tri    - Connectivity List
%    p      - Nodal Points
%
% Outputs:
%    plist  - NaN seperated polygon point list
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
% Get boundary edge segments of triangulation in the form xb (2xm) and
% yb (2xm). Each row corresponds the end points of each edge segment.
%---------------------------------------------------------------------
trep = triangulation(tri, p(:,1), p(:,2)); % Re-triangulate
fb = freeBoundary(trep); % Get boundary edge node pairs list

x = [p(fb(:,1),1) p(fb(:,2),1)]'; y = [p(fb(:,1),2) p(fb(:,2),2)]';

[ly,lx] = size(x);

%------------------------------------------------------------------------------
% Order boundary points
%------------------------------------------------------------------------------
% Initialize cell
plist = cell(numel(x),1);

% Store edge 1
plist{1} = [x(1,1) y(1,1)];
plist{2} = [x(2,1) y(2,1)];

ptr = 2; % Pointer to current end point

% Mark recorded points as nan's
x(1,1) = nan; y(1,1) = nan;
x(2,1) = nan; y(2,1) = nan;

i = 1; % Keep track of the start of each segment

Looping = 1;

while Looping
    
    % Find the next connecting edge with same end point
    ind = find( plist{ptr}(1) == x & plist{ptr}(2) == y ,1,'first');
    
    if isempty(ind) % If we're empty, are we done or moving on to the next boundary?
        
        % Close polygon
        ptr = ptr+1;
        plist{ptr} = plist{i};
                
        % Check to see if we're done
        if all(isnan(x(:))); Looping = 0; continue; end
        
        % Add a NaN to seperate
        ptr = ptr+1;
        plist{ptr} = [nan nan];
        
        ind = find(~isnan(x(1,:)),1,'first'); % Find a point to start
        
        ptr = ptr+1;
        plist{ptr} = [x(1,ind) y(1,ind)];
        
        % Store index to beginning of next point list
        i = ptr;
        
        ptr = ptr+1;
        plist{ptr} = [x(2,ind) y(2,ind)];
        
        % Mark recorded points as nan's
        x(1,ind) = nan; y(1,ind) = nan;
        x(2,ind) = nan; y(2,ind) = nan;
        
        continue; % Continue with loop
    end
    
    % convert linear index to subscript
    [r,c] = ind2sub([ly,lx],ind);
    
    if r == 1 % if r == 1, we want to grab the opposite end point
        
        ptr = ptr+1;
        plist{ptr} = [x(2,c) y(2,c)];
        
    else
        
        ptr = ptr+1;
        plist{ptr} = [x(1,c) y(1,c)];
        
    end
    
    % Mark recorded points as nan's
    x(1,c) = nan; y(1,c) = nan;
    x(2,c) = nan; y(2,c) = nan;
    
end

% Convert cell 2 matrix
plist = cell2mat(plist);

% Remove any potential double points
plist = unique(plist,'rows','stable');

% Close each boundary segment

% Create cell seperated segments
i1 = all(~isnan(plist),2); i2 = i1(:)';
idx = [strfind([~i2(1),i2],[0 1]); strfind([i2, ~i2(end)],[1 0])];
C = mat2cell(plist (i1,:),diff(idx)+1,size(plist,2));

% Loop over each segment and close boundary
for i = 1:length(C)
    
    C{i} = [C{i} ; C{i}(1,:); nan nan];
    
end

% Convert back to matrix
plist = cell2mat(C);

end