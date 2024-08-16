function MA_new = ConstructMainStream_v2(MA)
%% ========================================================================
% Construct main stream
%==========================================================================
Size = MA.Size;
ID = MA.ID;
ID1 = MA.ID;
Conn = MA.Conn;
Conn1 = MA.Conn;
XY = MA.XY;
id_channels = [];
ICONN = [];

ID_new = ID;
Conn_new = [];
XY_new = XY;
max_length = zeros(length(ID_new),1);
i = 0;
while 1
    if i >= size(Conn,1)
        break;
    end
    i = i + 1;
    
    % skip if both ends are not a tip
    if nnz(Conn == Conn(i,1)) ~= 1 && nnz(Conn == Conn(i,2)) ~= 1 
        continue;
    end
    
    num_channels = 0;
    id_channels = [];
    
    % Make the array "Conn" stores id in left-to-right order
    if nnz(Conn == Conn(i,2)) == 1
        Conn(i,:) = Conn(i,[2 1]);
%         XY{i} = flip(XY{i});
    end
    
    if nnz(Conn == Conn(i,2)) == 1 % skip if tip-to-tip
        continue;
    end
    
    num_channels = num_channels + 1;
    id_channels{num_channels} = Conn(i,:);
    
    Conn(i,:) = nan;
%     ID(i) = nan;
    k = 0;
    while 1
        k = k + 1;
        if k > length(id_channels)
            break;
        end
        ic = id_channels{k};
        
        j1 = find(Conn(:,1) == ic(end));
        j2 = find(Conn(:,2) == ic(end));
        Conn(j2,:) = Conn(j2,[2 1]);
%         for jj = 1 : length(j2)
%             XY{j2(jj)} = flip(XY{j2(jj)});
%         end
        j = vertcat(j1,j2);
        
        for jj = 1 : length(j)
            jjj = j(jj);
            if Conn(jjj,2) > 0
                id_channels{end+1} = [ic,Conn(jjj,2)];
            else
                id_channels{end+1} = [ic,Conn(jjj,2),nan];
            end
        end
        Conn(j,:) = nan;
%         ID(j) = nan;
    end
    
    % Remove channel connections which are entirely included in others
    id_channels = flip(id_channels);
    iremove = [];
    for j = 1 : length(id_channels)
        for k = j+1 : length(id_channels)
            if all(ismember(id_channels{k},id_channels{j}))
                id_channels{k} = [];
                iremove = [iremove,k];
            end
        end
    end
    id_channels(iremove) = [];

    % Compute length of channels
    LENGTH = zeros(length(id_channels),1);
    ICONN = cell(length(id_channels),1);
    for j = 1 : length(id_channels)
        id = id_channels{j};
        id(isnan(id)) = [];
        ICONN{j} = [];
        for k = 1 : (length(id)-1)
            kk = find(ismember(Conn1,id([k k+1]),'rows'));
            if isempty(kk)
                kk = -find(ismember(Conn1,id([k+1 k]),'rows'));
            end
            ICONN{j} = [ICONN{j}, kk];
            id1 = ID1{abs(kk)};
            [I,J] = ind2sub(Size,id1);
            dist = sum(sqrt( (I(1:end-1) - I(2:end)).^2 + (J(1:end-1) - J(2:end)).^2));
            LENGTH(j) = LENGTH(j) + dist;
        end
    end
    
    [~,id_main] = max(LENGTH);
    XY_add = [];
    id_new = [];
    ICONN_main = ICONN{id_main};
    for j = 1 : length(ICONN_main)
        if ICONN_main(j) > 0
            id_add = ID1{ICONN_main(j)};
            XY_add = vertcat(XY_add,XY{ICONN_main(j)});
        else
            id_add = flip(ID1{-ICONN_main(j)});
            XY_add = vertcat(XY_add,flip(XY{-ICONN_main(j)}));
        end
        id_new = vertcat(id_new,id_add);
        ID_new{abs(ICONN_main(j))} = [];
        
        
        
    end
    if ~isempty(id_new)
        ID_new{end+1} = id_new;
        XY_new{end+1} = XY_add;
        max_length(end+1) = max(LENGTH);
    end
    
%     ID_new = vertcat(ID{ICONN{id_main}});
    
%     ID{ICONN{id_main}(1)} = ID_new;
%     Conn1(ICONN{id_main}(1),:) = [ID_new(1), ID_new(end)];
    
%     Conn1(ICONN{id_main}(2:end),:) = 0;
    
    0;
    
    
end

id_remove = [];
for i = 1 : length(ID_new)
    ID_new{i} = unique(ID_new{i},'stable');
    XY_new{i} = unique(XY_new{i},'rows','stable');
    if isempty(ID_new{i})
        id_remove = [id_remove,i];
    end
end
ID_new(id_remove) = [];
XY_new(id_remove) = [];

for i = 1 : length(ID_new)
    Conn_new(i,:) = [ID_new{i}(1),ID_new{i}(end)];
end

% XY_new = cell(length(ID_new),1);
% for i = 1 : length(ID_new)
%     XY_new{i} = [X(ID_new{i}),Y(ID_new{i})];
% end

MA_new.nBranch = length(ID_new);
MA_new.ID = ID_new;
MA_new.Conn = Conn_new;
MA_new.XY = XY_new;


