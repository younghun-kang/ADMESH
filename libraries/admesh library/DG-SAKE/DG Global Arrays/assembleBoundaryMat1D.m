function [B,IL,IR,ZL,ZR] = assembleBoundaryMat1D(MESH,IL,IR,Z,b,h,ndof,npoints,ngauss,nelems)
% *************************************************************************
% Assembly/Compute the global advection matrices
%
% Author:
% Dustin West
% Assistant Engineer | Hazen and Sawyer
% dwest@hazenandsawyer.com | hazenandsawyer.com
% *************************************************************************

%------------------------------------------------------------------
% Initialize variables
%------------------------------------------------------------------
E       = MESH.edges;               % edge list
EA      = MESH.edgeAttachments;     % edge attachments
ET      = MESH.edgeType;            % edge attachments

nbp     = size(b.ii{1},2);          % # of 2D boundary points

IL.ix   = cell(nelems,1);           % Cell for 2D to 1D converging flow
IR.ix   = cell(nelems,1);           % Cell for 2D to 1D converging flow

% Also store left and right side slop values for flux
ZL      = cell(nelems,1);
ZR      = cell(nelems,1);

% Initialize global boundary flux matrix and sparse matrix indexing
[B,ib,jb] = deal(zeros(nelems*ndof*npoints,1));

% Initialize indexing vector for ib & jb
nId = ndof; id = 1:nId;

% Initialize counter for h
i = 1;

%------------------------------------------------------------------
% Compute some geometry
%-----------------------------------------------------------------

nelem2D = size(MESH.ConnectivityList,1);

% Assign element coordinates
xe = reshape(MESH.Points(MESH.ConnectivityList,1),nelem2D,3)';
ye = reshape(MESH.Points(MESH.ConnectivityList,2),nelem2D,3)';

% Compute edge lengths
L{1}    = sqrt((xe(3,:)-xe(2,:)).^2 + (ye(3,:)-ye(2,:)).^2);
L{2}    = sqrt((xe(1,:)-xe(3,:)).^2 + (ye(1,:)-ye(3,:)).^2);
L{3}    = sqrt((xe(2,:)-xe(1,:)).^2 + (ye(2,:)-ye(1,:)).^2);

% Compute edge normals
nx = zeros(3,nelem2D);                ny    = zeros(3,nelem2D);
nx(1,:) = (ye(3,:)-ye(2,:))./L{1}; ny(1,:) = (xe(2,:)-xe(3,:))./L{1};
nx(2,:) = (ye(1,:)-ye(3,:))./L{2}; ny(2,:) = (xe(3,:)-xe(1,:))./L{2};
nx(3,:) = (ye(2,:)-ye(1,:))./L{3}; ny(3,:) = (xe(1,:)-xe(2,:))./L{3};

clear xe ye L

%------------------------------------------------------------------
% Loop over each string in cascadeNetwork
%------------------------------------------------------------------
for k = 1:length(MESH.channelNetwork)
    
    % Assign global 1D node number for cascade
    ni = MESH.channelNetwork{k}(:,2);
    
    % Assign the 2D edges
    edges = [MESH.channelNetwork{k}(1:end-1,1) MESH.channelNetwork{k}(2:end,1)];
    
    % Number of elements in segment
    ne = size(edges,1);
    
    % Loop over each element in segment k
    for j = 1:ne
        
        % Place "left" element boundary matrix in the correct global slot
        B(id) = b.i(:,1)*(2./h(i));
        
        % Populate sparse matrix indexing
        ib(id) = (1:ndof)' + ndof*(i-1);
        jb(id) = ones(ndof,1)*ni(j);
        
        % Increment boundary index id
        id = id + nId;
        
        % Place "right" element boundary matrix in the correct global slot
        % Based on boundary condition
        B(id) = b.i(:,2)*(2./h(i));
        
        % Populate sparse matrix indexing
        ib(id) = (1:ndof)' + ndof*(i-1);
        jb(id) = ones(ndof,1)*ni(j+1);
        
        % Increment boundary index id
        id = id + nId;
        
        %----------------------------------------------------------------
        % Generate 2D to 1D index map for element
        %----------------------------------------------------------------
        
        % Find the row (edge number) in "E" that matches with the current
        % 1D edge we are looking at
        loc = find(...
            (edges(j,1) == E(:,1) | edges(j,1) == E(:,2)) & ...
            (edges(j,2) == E(:,1) | edges(j,2) == E(:,2)));
        
        % Retrieve "left" element attached to edge
        iL = EA{loc}(1);
        
        % Determine local edge number of edge (E(LOC)) from within element
        % iL
        k1 = find(MESH.ConnectivityList(iL,:)==E(loc,1));
        k2 = find(MESH.ConnectivityList(iL,:)==E(loc,2));
        kL = abs(k1-k2) + mod(max(k1,k2),3);
        
        % Retrieve "right" element attached to edge
        iR = EA{loc}(2);
        
        % Determine local edge number of edge (E(LOC)) from within element
        % iR
        k1 = find(MESH.ConnectivityList(iR,:)==E(loc,1));
        k2 = find(MESH.ConnectivityList(iR,:)==E(loc,2));
        kR = abs(k1-k2) + mod(max(k1,k2),3);
        
        % Compute side slope
        [zL,zR] = sideSlope(MESH,iL,iR,nx,ny,kL,kR,E(loc,:));
        
        % Store side slope
        ZL{i}    = ones(ngauss,1)*zL;
        ZR{i}    = ones(ngauss,1)*zR;
        
        % Check the edge type. We only want to store indeices where we have
        % a converging edge (thus the need to generate vector indexing for
        % lateral inflow). If the edte type is not == -2 then the edge is a
        % routing the flow so we don't need lateral inflow on these 1D
        % elements.
        if ET(loc) ~= -2
            
            % Enter zeros in place of left and right indexing values
            IL.ix{i} = zeros(ngauss,1);
            IR.ix{i} = zeros(ngauss,1);
            
            i = i + 1; % Move to next element
            continue;
        end
        
        % Assemble 2D to 1D indexing map-----------------------------------
        
        % Determine the 2D edge to global integration point index
        aL = ((iL-1)*3*nbp + (kL-1)*nbp + 1);
        
        % Generate index vector
        leftVec = (aL:(aL+nbp-1))';
        
        % Check if index needs to be fliped (order from highest to lowest
        if Z.edge2(leftVec(1)) < Z.edge2(leftVec(end))
            leftVec = flipud(leftVec);
        end
        
        % Find the location in IL.ii where indices leftVec are
        [~,IL.ix{i}] = ismember(leftVec,IL.ii);
        
        % Note, IL.ix, is going to be a "left" index vector to our 2D Fhat
        
        % Determine the edge to global integration point index
        aR = ((iR-1)*3*nbp + (kR-1)*nbp + 1);
        
        % Generate index vector
        rightVec = (aR:(aR+nbp-1))';
        
        % Check if index needs to be fliped
        if Z.edge2(rightVec(1)) < Z.edge2(rightVec(end))
            rightVec = flipud(rightVec);
        end
        
        % Find the location in IR.ii where indices rightVec are
        [~,IR.ix{i}] = ismember(rightVec,IR.ii);
        %------------------------------------------------------------------
        
        % Move to next element
        i = i + 1;
        
    end
    
end

% Convert cell to mat
IR.ix   = vertcat(IR.ix{:});
IL.ix   = vertcat(IL.ix{:});
ZL      = vertcat(ZL{:});
ZR      = vertcat(ZR{:});

% Remove zeros
B = B(ib > 0);
jb  = jb(ib > 0);
ib  = ib(ib > 0);

% Construct sparse boundary matrix
B   = sparse(ib,jb,B,nelems*ndof,npoints);

%==========================================================================
% SUBROUTINE
%==========================================================================

    function [zL,zR] = sideSlope(MESH,iL,iR,nx,ny,kL,kR,E)
        %
        % Author:
        % Dustin West
        % Assistant Engineer | Hazen and Sawyer
        % dwest@hazenandsawyer.com | hazenandsawyer.com
        
        % Compute edge length
        EL = sqrt( diff(MESH.Points(E,1))^2 + diff(MESH.Points(E,2))^2);
        
        % Compute mid-point of edge
        M.x = sum(MESH.Points(E,1))/2;
        M.y = sum(MESH.Points(E,2))/2;
        M.z = sum(MESH.Points(E,3))/2;
        
        % Compute "left" side slope----------------------------------------
        
        % Get coordinates of "left" element
        x = MESH.Points(MESH.ConnectivityList(iL,:),1);
        y = MESH.Points(MESH.ConnectivityList(iL,:),2);
        v = MESH.Points(MESH.ConnectivityList(iL,:),3);
        
        % Construct interpolant
        F = scatteredInterpolant(x,y,v,'linear','linear');
        
        % Project point into the element by a quarter of the edge length
        Q.x = M.x - (1/4)*EL*nx(kL,iL);
        Q.y = M.y - (1/4)*EL*ny(kL,iL);
        Q.z = F(Q.x,Q.y);
        
        % Compute distance between points
        d = sqrt( (M.x-Q.x)^2 + (M.y-Q.y)^2 );
        
        % Compute difference in elevation
        dz = abs(M.z - Q.z);
        
        % Compute angle assuming height of right sided triangle is 1
        theta = asin(dz/d);
        
        % Compute left side slope
        zL = 1/tan(theta);
        
        % Compute "right" side slope----------------------------------------
        
        % Get coordinates of element
        x = MESH.Points(MESH.ConnectivityList(iR,:),1);
        y = MESH.Points(MESH.ConnectivityList(iR,:),2);
        v = MESH.Points(MESH.ConnectivityList(iR,:),3);
        
        % Construct the interpolant
        F = scatteredInterpolant(x,y,v,'linear','linear');
        
        % Project point into the element by a quarter of the edge length
        Q.x = M.x - (1/4)*EL*nx(kR,iR);
        Q.y = M.y - (1/4)*EL*ny(kR,iR);
        Q.z = F(Q.x,Q.y);
        
        % Compute distance between points
        d = sqrt( (M.x-Q.x)^2 + (M.y-Q.y)^2 );
        
        % Compute difference in elevation
        dz = abs(M.z - Q.z);
        
        % Compute angle assuming height of right sided triangle is 1
        theta = asin(dz/d);
        
        % Compute left side slope
        zR = 1/tan(theta);
        
        %hold on
        %plot3([mx p(1)],[my p(2)],[mz p(3)],'r-o')
        
    end

end