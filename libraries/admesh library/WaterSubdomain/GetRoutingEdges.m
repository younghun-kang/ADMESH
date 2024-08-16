function GetRoutingEdges(app)
%--------------------------------------------------------------------------
% Load MESH from app and make a copy
%--------------------------------------------------------------------------
MESH0 = app.MESH;
MESH = MESH0;

%--------------------------------------------------------------------------
% Check if elements with zero-slope exist and modify elevation if so
%--------------------------------------------------------------------------
MESH = CheckElementSlope(MESH);

%--------------------------------------------------------------------------
% Flip z-values as functions developed for MATLAB version DG-SAKE use
% positive z value for topography
%--------------------------------------------------------------------------
MESH.Points(:,3) = -MESH.Points(:,3);

%--------------------------------------------------------------------------
% Change field name. ADMESH functions use 'Constraints.nodes' whereas
% DG-SAKE functions use 'Constraints.nodeStr'.
%--------------------------------------------------------------------------
[MESH.Constraints.nodes] = MESH.Constraints.nodeStr;
MESH.Constraints = rmfield(MESH.Constraints,'nodeStr');

%--------------------------------------------------------------------------
% Get routing edges
%--------------------------------------------------------------------------
MESH = fillSinks(MESH);                     % Fill sinks in mesh
MESH = assignEdgeType(MESH,0);              % Classify element edge type
MESH = routeDisjointedChannels(MESH);       % Route disconnected channels
MESH = assembleChannelNetwork(MESH);        % Assemble channel network
MESH = mesh_neighbors(MESH);

%--------------------------------------------------------------------------
% Return z-values to original sign
%--------------------------------------------------------------------------
MESH.Points(:,3) = -MESH.Points(:,3);

%--------------------------------------------------------------------------
% Construct basis functions
%--------------------------------------------------------------------------
p1 = 1;
p2 = 1;
[phi,psi,w] = DG_Element_Functions(p2,p1,MESH);

%--------------------------------------------------------------------------
% Construct DG element matrices
%--------------------------------------------------------------------------
fprintf('Assembling DG Element Arrays...\n\n');
[a,b,c] = DG_Element_Matrices(phi,w);

%--------------------------------------------------------------------------
% Assemble DG global arrays
%--------------------------------------------------------------------------
fprintf('Assembling DG Global Arrays...\n\n');
[A,B,C,IL,IR,BC,X,Y,Z,PHI,N,ZL,ZR] = DG_Global_Assembly(a,b,c,phi,psi,MESH);

%--------------------------------------------------------------------------
% Compute Kinematic Wave Properties
%--------------------------------------------------------------------------
fprintf('Constructing DG flux functions...\n\n');
n.elem2 = zeros(numel(X.elem2),1);
n.edge2 = zeros(numel(X.edge2),1);
n.elem1 = zeros(numel(X.elem1),1);
[F,G,ck,So] = GetKinematicWaveFlux(N,X,Y,ZL,ZR,MESH,n);

%--------------------------------------------------------------------------
% Load constraint data
%--------------------------------------------------------------------------
I = find([MESH.Constraints.num] == 18 | [MESH.Constraints.num] == 6);
nodes = vertcat(MESH.Constraints(I).nodes);
try
    data = vertcat(MESH.Constraints(I).data);
    data = data(:,2); % Take top width
catch
    data = ones(size(nodes,1),1);
end

%--------------------------------------------------------------------------
% Remove duplicated values (for degree of freedom in DG-SAKE)
%--------------------------------------------------------------------------
it = F.elem1.it(1:3:end);
ir = F.elem1.ir(1:3:end);

ZLR = ZL+ZR;
ZLR = ZLR(1:3:end);

%--------------------------------------------------------------------------
% Find nodes for channel network
%--------------------------------------------------------------------------
ncedno = cellfun(@(x) [x(1:end-1,1), x(2:end,1)],MESH.channelNetwork,'UniformOutput',0);
ncedno = vertcat(ncedno{:});

%--------------------------------------------------------------------------
% Check if a node is on a rectangular or triangular edge
%--------------------------------------------------------------------------
ncedno1 = unique(ncedno);
ir1 = zeros(size(ncedno1));
for i = 1 : length(ncedno1)
     ii = ncedno1(i);
     [I,J] = find(ncedno == ii);
     if any(ir(I))
         ir1(i) = 1;
     end
end

%--------------------------------------------------------------------------
% Set channel width
%--------------------------------------------------------------------------
channelwidth = zeros(length(ncedno1),2);
for i = 1 : length(ncedno1)
    ii = ncedno1(i);
    id = find(ncedno(:,1) == ii | ncedno(:,2) == ii);

    %----------------------------------------------------------------------
    % Use top-width taken from constraint data for rectangular channel, and
    % use averaged width from ZL and ZR for triangular channels.
    %----------------------------------------------------------------------
    if ir1(i)
        j = find(nodes == ii);
        channelwidth(i,:) = mean(data(j));
        if length(unique(data(j))) > 1
            warning(['Multiple values are found for a single node in a rectangular channel. ' ...
                'Averaged value is assigned.']);
        end
    else
        channelwidth(i,1) = mean(ZLR(id));
    end    
end

%--------------------------------------------------------------------------
% Reorganize channel width
%--------------------------------------------------------------------------
ChannelData = cellfun(@(x) x(:,1),MESH.channelNetwork,'UniformOutput',0);
for i = 1 : length(ChannelData)
    cd = zeros(length(ChannelData{i,1}),2); 
    for j = 1 : length(ChannelData{i,1})
        I = (ncedno1 == ChannelData{i,1}(j));
        cd(j,:) = channelwidth(I,:);
    end
    ChannelData{i,2} = cd;
end

%--------------------------------------------------------------------------
% Remove triangular channels that overlap rectangular channels
%--------------------------------------------------------------------------
for i = 1 : size(ChannelData,1)
    j = find(ChannelData{i,2}(:,2) > 0,1,'first');
    if isempty(j)
        continue;
    end
    I = 1:j;
    ChannelData{i,1} = ChannelData{i,1}(I,:);
    ChannelData{i,2} = ChannelData{i,2}(I,:);
end

%--------------------------------------------------------------------------
% Remove empty or single-point triangular channels
%--------------------------------------------------------------------------
I = cellfun(@(x) length(x),ChannelData(:,1));
ChannelData = ChannelData(I > 1,:);

%--------------------------------------------------------------------------
% Return field name to ADMESH style
%--------------------------------------------------------------------------
[MESH.Constraints.nodeStr] = MESH.Constraints.nodes;
MESH.Constraints = rmfield(MESH.Constraints,'nodes');

%--------------------------------------------------------------------------
% Add routing edges as IBTYP=6
%--------------------------------------------------------------------------
N = length(MESH.Constraints);
for i = 1 : size(ChannelData,1)
    MESH.Constraints(N+i).num = 6;
    MESH.Constraints(N+i).nodeStr = ChannelData{i};
    MESH.Constraints(N+i).data = [ones(length(ChannelData{i}),1), ChannelData{i,2}];
end

app.MESH = MESH;

PlotMesh(app,.1);





