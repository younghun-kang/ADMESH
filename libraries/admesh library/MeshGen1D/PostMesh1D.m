function [Mesh1D,Constraints] = PostMesh1D(Mesh1D,Constraints,hmin)
Mesh1Dorigin = Mesh1D; % Keep original Mesh1D for debugging
Mesh1Dpost = cell(0);  % Save step-by-step Mesh1D for debugging

%% ========================================================================
% Transfer shoreline constraints to internal boundaries if it is straightline
%==========================================================================
IC3 = find(abs([Constraints.num]) == 19);
if ~isempty(IC3)
Mesh1D_c3 = Mesh1D(IC3);
IC3toC2 = [];
for i = 1 : length(Mesh1D_c3)
    if isequal([Mesh1D_c3(i).X(1),Mesh1D_c3(i).Y(1)],[Mesh1D_c3(i).X(end),Mesh1D_c3(i).Y(end)])...
        && length(unique([Mesh1D_c3(i).X,Mesh1D_c3(i).Y],'rows')) == 2
        IC3toC2 = [IC3toC2;IC3(i)];
    end
end

for i = 1 : length(IC3toC2)
    j = IC3toC2(i);
    Constraints(j).num = sign(Constraints(j).num)*17;
end
end
Mesh1Dpost{end+1} = Mesh1D;

%% ========================================================================
% Merge 1D point clusters which have neighboring fixed points within hmin/2
%==========================================================================
Points1D = [vertcat(Mesh1D.X), vertcat(Mesh1D.Y)];

% Find close points
[id,dist] = knnsearch(Points1D,Points1D,'k',10);

% Check the number of points closer than hmin/4
idClose = dist < hmin/4;
numClose = sum(idClose,2);

% Find cluster IDs
I = find(numClose > 1);
id1 = cell(length(I),1);
for j = 1 : length(I)
    i = I(j);
    id1{j} = sort(id(i,1:numClose(i)))';
end
id1 = id1(:);

% Remove duplicate cluster
for i = 1 : length(id1)
    for j = i+1 : length(id1)
        if all(ismember(id1{j},id1{i}))
            id1{j} = [];
        end
        if all(ismember(id1{i},id1{j}))
            id1{i} = [];
        end
    end
end
I = cellfun(@(x) ~isempty(x),id1); 
id1 = id1(I);

% Merger clusters if overlapping clusters exist
for i = 1 : length(id1)
    for j = i+1 : length(id1)
        if any(ismember(id1{j},id1{i}))
            temp = vertcat(id1{j},id1{j});
            [~,dist1] = knnsearch(Points1D(temp,:),Points1D(temp,:),'k',length(temp));
            if max(dist1(:)) > hmin/2
                warning('There are overlapping clusters. These are merged with diameter %f.',max(dist1(:)));
            end
            id1{i} = unique(vertcat(id1{i},id1{j}));
            id1{j} = [];
        end
    end
end
I = cellfun(@(x) ~isempty(x),id1); 
id1 = id1(I);

% Check cluster diameters
clusterD = zeros(length(id1),1);
for i = 1 : length(id1)
    [~,temp] = knnsearch(Points1D(id1{i},:),Points1D(id1{i},:),'k',length(id1{i}));
    clusterD(i) = max(temp(:));
end
if max(clusterD/hmin) > 0.5
    warning('Some clusters have diameter greater than 0.5*hmin (%.2f).',max(clusterD/hmin));
end

mfp = [];
iMap = [];
for i = 1 : length(id1)
    mfp(i,:) = mean(Points1D(id1{i},:));
    iMap = [iMap; i*ones(length(id1{i}),1)];
end
mfp(:,3) = 0;

iFixedMerge = vertcat(id1{:});
for i = 1 : length(Mesh1D)
    x = Mesh1D(i).X;
    y = Mesh1D(i).Y;
    K = Mesh1D(i).K;

    [I,J] = ismember([x,y],Points1D(iFixedMerge,:),'rows');

    if any(I)
        I = find(I);
        J = nonzeros(J);
        J = iMap(J);

        x(I) = mfp(J,1);
        y(I) = mfp(J,2);
        K(I) = mfp(J,3);

        k1 = find(diff(I) == 1);
        k2 = I(k1);
        k = k2(J(k1) == J(k1+1));
        x(k) = [];
        y(k) = [];
        K(k) = [];
    end
    Mesh1D(i).X = x;
    Mesh1D(i).Y = y;
    Mesh1D(i).K = K;
end

I = cellfun(@(x) length(x),{Mesh1D.X}) > 1;
Mesh1D = Mesh1D(I);
Constraints = Constraints(I);

Mesh1Dpost{end+1} = Mesh1D;

%% ========================================================================
% Remove 1D elements nodes with length smaller than hmin/2
%==========================================================================

fixedPoints = cell(length(Mesh1D),1);
for i = 1 : length(Mesh1D)
    x = Mesh1D(i).X;
    y = Mesh1D(i).Y;
    fixedPoints{i} = [x([1 end]), y([1 end])];
end
fixedPoints = vertcat(fixedPoints{:});

for i = 1 : length(Mesh1D)
    x = Mesh1D(i).X;
    y = Mesh1D(i).Y;
    k = Mesh1D(i).K;
    while 1 % Remove 1D element shorter than hmin/2 except fixed points
        d = sqrt(diff(x).^2 + diff(y).^2);
        I = find(d < hmin/2);

        iFixed = find(ismember([x,y],fixedPoints,'rows'));
        I = setdiff(I,iFixed);
        if isempty(I)
            break;
        end
        [~,J] = min(d(I));
        I = I(J);
        x(I(1)) = [];
        y(I(1)) = [];
        k(I(1)) = [];
    end
    Mesh1D(i).X = x;
    Mesh1D(i).Y = y;
    Mesh1D(i).K = k;
end
Mesh1Dpost{end+1} = Mesh1D;

%% ========================================================================
% Remove 1D constraints that completely included in shoreline constraints
%==========================================================================
IC3 = find(abs([Constraints.num]) == 19);
IC3 = IC3(:);
if isempty(IC3)
    return;
end
C3XY = [vertcat(Mesh1D(IC3).X),vertcat(Mesh1D(IC3).Y)];

for i = 1 : length(Mesh1D)
    if abs(Constraints(i).num) == 19
        continue;
    end
    x = Mesh1D(i).X;
    y = Mesh1D(i).Y;
    K = Mesh1D(i).K;

    [~,dist] = knnsearch(C3XY,[x,y]);

    if max(dist) < hmin/2
        Mesh1D(i).X = [];
        Mesh1D(i).Y = [];
        Mesh1D(i).K = [];
    end
end

I = cellfun(@(x) length(x),{Mesh1D.X}) > 1;
Mesh1D = Mesh1D(I);
Constraints = Constraints(I);
Mesh1Dpost{end+1} = Mesh1D;

%% ========================================================================
% Find open-channels/internal boundaries too close to shoreline using VDT
%==========================================================================
IC3 = find(abs([Constraints.num]) == 19);
IC3 = IC3(:);
C3X = {Mesh1D(IC3).X};
C3Y = {Mesh1D(IC3).Y};
C3XY = cellfun(@(x,y) [x(:),y(:)],C3X,C3Y,'UniformOutput',0);

dx = 1; % resolution of background grid (convert meter to degree)
temp2 = 0;
for i = 1 : length(Mesh1D)
    if abs(Constraints(i).num) == 19
        continue;
    end
    x = Mesh1D(i).X;
    y = Mesh1D(i).Y;
    K = Mesh1D(i).K;

    [~,~,dist] = Compute8SSED_v3(C3XY,x,y,dx);

    while 1
        Mesh1Dcopy = Mesh1D;
        if dist(1) > 0 && dist(1) < hmin/10
            Mesh1D(i).X(1) = [];
            Mesh1D(i).Y(1) = [];
            Mesh1D(i).K(1) = [];
            dist(1) = [];
        end
        if dist(end) > 0 && dist(end) < hmin/10
            Mesh1D(i).X(end) = [];
            Mesh1D(i).Y(end) = [];
            Mesh1D(i).K(end) = [];
            dist(end) = [];
        end
        if isequal(Mesh1D,Mesh1Dcopy)
            break;
        else
            temp2 = temp2+1;
        end
    end
end

I = cellfun(@(x) length(x),{Mesh1D.X}) > 1;
Mesh1D = Mesh1D(I);
Constraints = Constraints(I);
Mesh1Dpost{end+1} = Mesh1D;

%% ========================================================================
% Find open-channels/internal boundaries crossing shorelines
%==========================================================================
IC3 = find(abs([Constraints.num]) == 19);
IC3 = IC3(:);
C3X = {Mesh1D(IC3).X};
C3Y = {Mesh1D(IC3).Y};
C3XY = cellfun(@(x,y) [x(:),y(:); nan, nan],C3X,C3Y,'UniformOutput',0);
C3XY = vertcat(C3XY{:});

temp2 = 0;
for i = 1 : length(Mesh1D)
    if abs(Constraints(i).num) == 19
        continue;
    end

    x = Mesh1D(i).X;
    y = Mesh1D(i).Y;

    [xi,yi] = polyxpoly(x,y,C3XY(:,1),C3XY(:,2));
    [~,dist] = knnsearch([x,y],[xi,yi]);

    if dist > 0
        % Remove this constraint as it's not clear which side should be removed at this point
        temp2 = temp2 + 1;
        Mesh1D(i).X = [];
        Mesh1D(i).Y = [];
        Mesh1D(i).K = [];
    end
end

I = cellfun(@(x) length(x),{Mesh1D.X}) > 1;
Mesh1D = Mesh1D(I);
Constraints = Constraints(I);

Mesh1Dpost{end+1} = Mesh1D;
