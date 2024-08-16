function [pq,tq] = tri2quad(p,t)

%-------------------------------------------------------------------------
% Get triangulation representation
%-------------------------------------------------------------------------
trep = triangulation(t,p(:,[1 2]));

%-------------------------------------------------------------------------
% Construct a set of edges that join the circumcenters of neighboring
% triangles; the additional logic constructs a unique set of such edges
%-------------------------------------------------------------------------
neigh = neighbors(trep);            % Get neighbor list
ic = incenter(trep);                % Find points in center
ctx = ic(:,1); cty = ic(:,2);       % Distribute to ctx & cty
numt = size(trep,1); T = (1:numt)'; % Number of elements

idx1 = T < neigh(:,1);
idx2 = T < neigh(:,2);
idx3 = T < neigh(:,3);

neigh = [T(idx1) neigh(idx1,1); T(idx2) neigh(idx2,2); T(idx3) neigh(idx3,3)];

%-------------------------------------------------------------------------
% Remove all connection edges that span constraints in our mesh
%-------------------------------------------------------------------------
% ci = find([MESH.Constraints.num] == 18);
% 
% for k = 1:numel(ci)
%     
%     c = MESH.Constraints(ci(k)).nodeStr;
% 
%     ti = sort(cell2mat(edgeAttachments(trep,c(1:end-1),c(2:end))),2);
%     
%     % Remove connections
%     neigh = setdiff(neigh',ti,'rows','stable')';
%     
% end

%-------------------------------------------------------------------------
% For each point (in center) we have up to 3 possibilities for merging
% elements. Pick the one that will give the highest element quality. For
% elements with no neighbors, split element into a quadrilateral.
%-------------------------------------------------------------------------
neigh = sort(neigh,2)'; % Sort connecting edges matrix
qt = [1 2 3 4];         % Simple connectivity used for quads
q = cell(1,3);          % Cell for holding up to three quad elements
Quad = cell(numt,1);    % Cell used for storing quadrilateral elements
Q = zeros(1,3);         % Initialize element quality vector

% Loop over each element
for k = 1:numt
    
    if mod(k,150) == 0
        uiStatusBar(k/numt)
    end
    
    % Assign nodes for element k
    xe = p(trep(k,:),1);
    ye = p(trep(k,:),2);
    
    % Find edge connections that include element k
    [r,c] = find(k == neigh);
    
    % No matches found, split triangular element
    if isempty(c);
                
        % Compute midpoint of each edge
        mx = (xe([1 2 3]) + xe([2 3 1]))./2;
        my = (ye([1 2 3]) + ye([2 3 1]))./2;
                
        % Create quad elements
        quadx = [...
            [ctx(k) mx(3) xe(1) mx(1)]';...
            [ctx(k) mx(1) xe(2) mx(2)]';...
            [ctx(k) mx(2) xe(3) mx(3)]'];
        
        quady = [...
            [cty(k) my(3) ye(1) my(1)]';...
            [cty(k) my(1) ye(2) my(2)]';...
            [cty(k) my(2) ye(3) my(3)]'];
        
        Quad{k} = [quadx quady];
        
        continue;
    end
    
    % Initialize element quality vector
    Q(:) = 0;
    
    % For each neighboring element, create quad & test quality
    for j = 1:numel(c)

        % Neighbor j
        if r(j) == 1
            nj = neigh(2,c(j));
        else
            nj = neigh(1,c(j));
        end
        
        % Coordinates to neighbor j
        xj = p(trep(nj,:),1);
        yj = p(trep(nj,:),2);
        
        % Find the points to join
        lie = ismember([xe,ye],[xj,yj],'rows');
        lij = ismember([xj,yj],[xe,ye],'rows');
        
        % Create quadrilateral
        q{j} = [ ...
            [xe(~lie)               ,ye(~lie)];...
            [xe(find(lie,1,'first')), ye(find(lie,1,'first'))];...
            [xj(~lij)               ,yj(~lij)];...
            [xe(find(lie,1,'last')), ye(find(lie,1,'last'))]];
        
        % Test quality
        [~,~,Q(j)] = MeshQuality(q{j},qt,0,'Quad');
                
    end
    
    % Keep connection with maximum element quality & remove other connections
    [~,maxq] = max(Q);                        % Find max Q
    neigh(:,c(maxq ~= (1:numel(c)))) = nan;   % Remove other connections to element k
    
    % Find other connections in new neighboring element
    [~,ind] = find(neigh(2,c(maxq)) == neigh);
    neigh(:,ind(ind~=c(maxq))) = nan; % Remove connections
    
    % Combine neighbors into a Quad element
    
    % Compute midpoint of each edge
    mx = (q{maxq}(qt(1,[1 2 3 4]),1) + q{maxq}(qt(1,[2 3 4 1]),1))./2;
    my = (q{maxq}(qt(1,[1 2 3 4]),2) + q{maxq}(qt(1,[2 3 4 1]),2))./2;
    
    % Compute centroid
    [~,cx,cy] = polycenter(q{maxq}(:,1),q{maxq}(:,2));
    
    % Store quads in cell array
    Quad{k} = [...
    [q{maxq}(1,1),q{maxq}(1,2); mx(1),my(1); cx,cy; mx(4),my(4); ];... % Quad 1
    [q{maxq}(2,1),q{maxq}(2,2); mx(2),my(2); cx,cy; mx(1),my(1); ];... % Quad 2
    [q{maxq}(3,1),q{maxq}(3,2); mx(3),my(3); cx,cy; mx(2),my(2); ];... % Quad 3
    [q{maxq}(4,1),q{maxq}(4,2); mx(4),my(4); cx,cy; mx(3),my(3); ];];  % Quad 4
        
end

%-------------------------------------------------------------------------
% Construct mesh structure for quad mesh
%-------------------------------------------------------------------------

% Convert structure to cell, distribute to x & y vecs
Quad = cell2mat(Quad);
quadx = reshape(Quad(:,1),4,size(Quad,1)/4)';
quady = reshape(Quad(:,2),4,size(Quad,1)/4)';

clear Quad

% Create a unique points list
pq =  unique([quadx(:),quady(:)],'rows','stable');

% Generate connectivity list
[~,lib1] = ismember([quadx(:,1) quady(:,1)],pq,'rows');
[~,lib2] = ismember([quadx(:,2) quady(:,2)],pq,'rows');
[~,lib3] = ismember([quadx(:,3) quady(:,3)],pq,'rows');
[~,lib4] = ismember([quadx(:,4) quady(:,4)],pq,'rows');

tq = [lib1 lib2 lib3 lib4];

clear lib1 lib2 lib3 lib4



end