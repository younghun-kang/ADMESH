SHP = shaperead('D:\academic\data\gis\walnut_gulch\SWRC\streamlines\streamlines.shp');
% SHP = struct2table(SHP);

%% ========================================================================
% Merge channels based on stream order field
%==========================================================================
ORDER = [6 5 4];
XY = cell(0);
data = cell(0);
for k = 1 : length(ORDER)
order = ORDER(k);

I = ([SHP.ORDER_] == order);
SHP1 = SHP(I,:);

%--------------------------------------------------------------------------
% Find joints
%--------------------------------------------------------------------------
X = {SHP1.X};
Y = {SHP1.Y};
X = cellfun(@(x) x(~isnan(x)),X,'UniformOutput',0);
Y = cellfun(@(x) x(~isnan(x)),Y,'UniformOutput',0);
X = cellfun(@(x) x(:),X,'UniformOutput',0);
Y = cellfun(@(x) x(:),Y,'UniformOutput',0);
X = vertcat(X{:});
Y = vertcat(Y{:});
[~,I] = unique([X(:),Y(:)],'rows');
X(I) = [];
Y(I) = [];
XYj = [X(:), Y(:)];

%--------------------------------------------------------------------------
% Find end points
%--------------------------------------------------------------------------
X = {SHP1.X};
Y = {SHP1.Y};
X = cellfun(@(x) x(~isnan(x)),X,'UniformOutput',0);
Y = cellfun(@(x) x(~isnan(x)),Y,'UniformOutput',0);
XY1 = cellfun(@(x,y) [x(1) y(1)],X,Y,'UniformOutput',0);
XY2 = cellfun(@(x,y) [x(end) y(end)],X,Y,'UniformOutput',0);
XY1 = vertcat(XY1{:});
XY2 = vertcat(XY2{:});

%--------------------------------------------------------------------------
% Find id of end points
%--------------------------------------------------------------------------
[~,I1] = ismember(XY1,XYj,'rows');
[~,I2] = ismember(XY2,XYj,'rows');
IN = [I1(:), I2(:)];
% I1 = nonzeros(I1);
% I2 = nonzeros(I2);

%--------------------------------------------------------------------------
% Find end of streamline
%--------------------------------------------------------------------------
[I,J] = find(I1 == 0);
IEs = num2cell(I);
JEs = num2cell(J);

%--------------------------------------------------------------------------
% Construct streamline
%--------------------------------------------------------------------------
for j = 1 : length(IEs)
    IE = IEs{j};
    JE = JEs{j};
    for i = 1 : size(IN,1)
        iE = IE(end);
        jE = JE(end);

        if jE == 1
            iN = IN(iE,2);
        elseif jE == 2
            iN = IN(iE,1);
        end

        if iN == 0
            continue;
        end

        [iE2,jE2] = find(IN == iN);
        I = iE2 == iE | jE2 == 2;
        jE2(I) = [];
        iE2(I) = [];

        if isempty(iE2)
            continue;
        end

        if length(iE2) > 1
            continue;
        end

        IE(end+1) = iE2;
        JE(end+1) = jE2;

        IN(IN==iN) = 0;

        if all(IN == 0)
            break;
        end
    end
    IEs{j} = IE;
    JEs{j} = JE;
end

if order == 6
    W = 33.2999;
    H = 0.7969;
elseif order == 5
    W = 19.6764;
    H = 0.5257;
elseif order == 4
    W = 9.2481;
    H = 0.5290;
end

for j = 1 : length(IEs)
    IE = IEs{j};
    X = {SHP1.X};
    Y = {SHP1.Y};
    X = cellfun(@(x) x(~isnan(x)),X,'UniformOutput',0);
    Y = cellfun(@(x) x(~isnan(x)),Y,'UniformOutput',0);
    X = cellfun(@(x) x(:),X,'UniformOutput',0);
    Y = cellfun(@(x) x(:),Y,'UniformOutput',0);
    X = vertcat(X{IE});
    Y = vertcat(Y{IE});
    [~,I] = unique([X(:), Y(:)],'rows','stable');
    XY{end+1} = [X(I), Y(I)];
    data{end+1} = repmat([H W W],length(I),1);
end

end

%% ========================================================================
% Plot result
%==========================================================================
figure; hold on;
for j = 1 : length(IEs)
    IE = IEs{j};
    X = {SHP1.X};
    Y = {SHP1.Y};
    X = cellfun(@(x) x(~isnan(x)),X,'UniformOutput',0);
    Y = cellfun(@(x) x(~isnan(x)),Y,'UniformOutput',0);
    X = cellfun(@(x) x(:),X,'UniformOutput',0);
    Y = cellfun(@(x) x(:),Y,'UniformOutput',0);
    X = vertcat(X{IE});
    Y = vertcat(Y{IE});

    plot(X,Y);
    plot(X(1),Y(1),'o');
end
axis equal;

figure; hold on;
for i = 1 : length(XY)
    x = XY{i}(:,1);
    y = XY{i}(:,2);
    z = data{i}(:,2);
    myPlot3(x,y,z);
end

%% ========================================================================
% Write ADMESH input file
%==========================================================================
SHP = shaperead('D:\academic\data\gis\walnut_gulch\SWRC\boundary\boundary.shp');

clear PTS;
PTS.Poly.x = SHP.X;
PTS.Poly.y = SHP.Y;

k = 0;
for i = 1 : length(XY)
    k = k + 1;
    PTS.Constraints(k).num = -18;
    PTS.Constraints(k).xy = XY{i};
    PTS.Constraints(k).type = 'line';
    PTS.Constraints(k).data = data{i};
    PTS.Constraints(k).Kappa = zeros(size(XY{i},1),1);
end

Settings.DummyConstraint = 1;

filename = ['test_files/GIS input/Walnut Gulch SWRC/WGEW_SWRC_stream',num2str(min(ORDER))];
save(filename,'PTS');

