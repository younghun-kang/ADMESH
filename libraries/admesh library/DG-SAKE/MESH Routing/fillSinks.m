function MESH = fillSinks(MESH)
%==========================================================================
% Remove sinks in MESH
%
% The fillsnks algorithm looks at each interior vertex (vertices not
% belonging to the boundary of the domain) and checks if the vertex is a
% sink. If the vertex is a sink the elevation of that vertex is elevated to
% the minimum of the neighboring vertices.
%
% One potential flaw in the current algorithm could lead to flats in 1D
% cascade segments. Additional work is needed in this function to
% accomadate this.
%
% Author:
% Dustin West
% Assistant Engineer | Hazen and Sawyer               
% dwest@hazenandsawyer.com | hazenandsawyer.com
%==========================================================================

display = 0;

fprintf('Filling sinks...\n\n');

% Get triangulation representation
trep = triangulation(MESH.ConnectivityList,MESH.Points);

% Get the edges corresponding to the boundaries
fb = freeBoundary(trep);

% Find indices corresponding to interior nodes only
int = setdiff( (1:size(MESH.Points,1))',fb);

% Find the elements attached to each vertex in "int"
vt = vertexAttachments(trep, int);

% Initialize original and 'filled' elevation vectors
z   = MESH.Points(:,3); % Original elevation
zf  = z;                % Filled elevation variable

% Initialize conditional flag
flag = true;

% Enter while loop
while flag == 1
    
    flag = 0; % Set flag to zero in case we have no sinks
    
    % Loop over each internal node
    for i = 1:length(int)
        
        % Find the vertices attached to vertex (va)
        ti = MESH.ConnectivityList(vt{i},:);
        va = ti(ti ~= int(i));
        
        % Compute the difference between node i and the surrounding va
        dz = (z(int(i)) - z(va));
        
        % Check for sink (true is z(in(i),3) is the lowest point)
        if all(dz < eps & abs(dz) > eps)
            %zf(int(i)) = mean(z(va)); % Set elevation to minimum
            %zf(int(i)) = min(z(va)) + min(z(va))*0.0005; % Set elevation to minimum
            %zf(int(i)) = min(z(va)) +  min(z(va))*0.0000009; % Set elevation to minimum
            zf(int(i)) = min(z(va)) +  min(z(va))*0.00006; % Set elevation to minimum
            %zf(int(i)) = min(z(va)) +  min(z(va))*0.0001; % Set elevation to minimum
            flag = 1;                % Reset flag
        end
        
    end
    
    % Assign new elevation
    z = zf;
    
end

% Store z in MESH data structure
MESH.Points(:,3) = z;


%==========================================================================
% For debugging
%==========================================================================
if display == 1
    lowest = false(size(z));
    
    for i = 1:length(int)
        
        % Find the vertices attached to vertex
        ti = MESH.ConnectivityList(vt{i},:);
        va = ti(ti ~= int(i));
        
        % Compute the difference
        dz = (z(int(i)) - z(va));
        
        % Check for sink (true is MESH.Points(in(i),3) is the lowest point)
        if all(dz < eps & abs(dz) > eps)
            lowest(i) = true;
        end
        
    end
    
    %--------------------------------------------------------------------------
    % Display results
    %--------------------------------------------------------------------------
    zmin = min(MESH.Points(:,3)); % Compute min elevation
    zmax = max(MESH.Points(:,3)); % Compute max elevation
    
    % Normalize and multiply by 256 (number of rows in colormap)
    C = round(256*(MESH.Points(:,3) - zmin) / ( zmax - zmin ));
    C(C>256) = 256;  C(C <= 0) = 1; % Cap
    
    % Get landcover colormap
    cmap = flipud(landcolor(256));
    
    %--------------------------------------------------------------------------
    % Plot MESH
    %--------------------------------------------------------------------------
    patch(...
        'faces',MESH.ConnectivityList,...
        'vertices',MESH.Points,...
        'facevertexcdata',cmap(C,:),...
        'facecolor','interp',...
        'edgecolor','k',...
        'edgealpha',.2);
    hold on
    x = MESH.Points(:,1);
    y = MESH.Points(:,2);
    %triplot(trep)
    plot3(x(fb'),y(fb'),z(fb'),'r')
    plot3(x(lowest),y(lowest),z(lowest),'m.','markersize',50)
    [ts,~,~] = triSlopeVector(MESH);
    
    % Compute centroid
    ic = (...
        MESH.Points(MESH.ConnectivityList(:,1),:) + ...
        MESH.Points(MESH.ConnectivityList(:,2),:) + ...
        MESH.Points(MESH.ConnectivityList(:,3),:))/3;
    
    % Compute elevation offset
    offset = ic(:,3)*.00001;
    
    % Display flow direction
    quiver3(ic(:,1),ic(:,2),ic(:,3)+offset,ts(:,1),ts(:,2),ts(:,3),.5,'color','k')
    
    daspect([1 1 1/23])
end

end