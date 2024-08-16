function ReorderOpenChannel(app)

%--------------------------------------------------------------------------
% Check if MESH exists
%--------------------------------------------------------------------------
if isempty(app.MESH)
    msg = 'No MESH data is found.';
    uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');
    return;
end

%--------------------------------------------------------------------------
% Get data from input
%--------------------------------------------------------------------------
MESH = app.MESH;
if any([MESH.Constraints.num] == 18)
    ix = find([MESH.Constraints.num] == 18);
    if ~isempty(app.CrossSection)
        CS = app.CrossSection;
        try
            data = cellfun(@(x) [x.H, x.TW, x.BW],CS,'UniformOutput',0);
        catch
            data = [];
        end

        if isempty(data)
            return;
        end

        for i = 1 : length(ix)
            if length(MESH.Constraints(ix(i)).nodeStr) ~= size(data{i},1)
                error;
            end
            MESH.Constraints(ix(i)).data = data{i};
            MESH.Constraints(ix(i)).num = 6;
        end
        % Nodal segment
    end

    if any(cellfun(@(x) size(x,2),{MESH.Constraints(ix).data}))
        for i = 1 : length(ix)
            if size(MESH.Constraints(ix(i)).data,2) == 3
                MESH.Constraints(ix(i)).num = 6;
            end
        end
    end
end

PTS = MESH.Points;
Constraints = MESH.Constraints;
Topo = -PTS(:,3);

%--------------------------------------------------------------------------
% Construct new internal boundaries for IBTYPE=6
%--------------------------------------------------------------------------
NewNodes = [];
NewData = [];
idType6 = find([Constraints.num]==6);
for k = 1 : length(idType6)
    i = idType6(k);

    Nodes = Constraints(i).nodeStr;
    Data = Constraints(i).data;
    NodesTopo = Topo(Nodes);

    %----------------------------------------------------------------------
    % Remove local maxima of channels with given threshold
    %----------------------------------------------------------------------
    while 1
        dNT = diff(NodesTopo);
        I = find(dNT >= 0 & dNT < inf) + 1;
        if isempty(I)
            break;
        end
        for j = 1 : length(I)
            if NodesTopo(I(j)-1) == 0
                NodesTopo(I(j)) = -0.00006;
            else
                NodesTopo(I(j)) = NodesTopo(I(j)-1)*(1 -sign(NodesTopo(I(j)-1))*0.00006);
            end
        end
    end
    Topo(Nodes) = NodesTopo;

    %----------------------------------------------------------------------
    % Split channels so that each channel has all downslopes
    %----------------------------------------------------------------------
    idLocMax = find(islocalmax(NodesTopo));
    idLocMin = find(islocalmin(NodesTopo));

    if NodesTopo(1) > NodesTopo(2)
        idLocMax = [1; idLocMax];
    else
        idLocMin = [1; idLocMin];
    end

    if NodesTopo(end) > NodesTopo(end-1)
        idLocMax = [idLocMax; length(NodesTopo)];
    else
        idLocMin = [idLocMin; length(NodesTopo)];
    end

    while 1
        if isempty(idLocMax) || isempty(idLocMin)
            break;
        end
        i = idLocMax(1);
        j = idLocMin(1);

        if i < j
            NewNodes{end+1} = Nodes(i:j);
            NewData{end+1} = Data(i:j,:);
            idLocMax(1) = [];
        else
            NewNodes{end+1} = Nodes(i:-1:j);
            NewData{end+1} = Data(i:-1:j,:);
            idLocMin(1) = [];
        end
    end
end

%--------------------------------------------------------------------------
% Remove original internal constraints with IBTYPE=6
%--------------------------------------------------------------------------
Constraints(idType6) = [];

%--------------------------------------------------------------------------
% Add new internal constraints
%--------------------------------------------------------------------------
n = length(Constraints);
for i = 1 : length(NewNodes)
    Constraints(n+i).num = 6;
    Constraints(n+i).nodeStr = NewNodes{i};
    Constraints(n+i).data = NewData{i};
end

%--------------------------------------------------------------------------
% Save new internal constraints
%--------------------------------------------------------------------------
MESH.Points(:,3) = -Topo;
app.MESH = MESH;
app.MESH.Constraints = Constraints;









