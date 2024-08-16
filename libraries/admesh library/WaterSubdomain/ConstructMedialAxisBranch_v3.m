function MA = ConstructMedialAxisBranch_v3(id_MA)
N = size(id_MA,1);
M = size(id_MA,2);

%--------------------------------------------------------------------------
% Identify joints based on the number of connected MA points
%--------------------------------------------------------------------------
[I,J] = find(id_MA);
id_joint = sparse(N,M);
for ii = 1 : length(I)
    i = I(ii);
    j = J(ii);
    k = N*(j-1) + i;
        
    if id_MA(k) == 0
        continue;
    end
    %         if ~ismember(k,idIN)
    %             continue;
    %         end
    %         k_nghb = k + [-M-1, -M, -M+1, -1, 1, M-1, M, M+1];
    k_nghb = k + [-N-1, -N, -N+1, -1, 1, N-1, N, N+1];
    temp = id_MA(k_nghb) == 1;
    
    if nnz(temp) > 2
        id_joint(k) = 1;
    end
end

id_joint = logical(id_joint);

%--------------------------------------------------------------------------
% Compute the number of connections
%--------------------------------------------------------------------------
id_MA2 = find(id_MA);
num_nghb_MA = zeros(length(id_MA2),1);
for i = 1 : length(id_MA2)
    k = id_MA2(i);
    k_nghb1 = k + [-N-1, -N, -N+1, -1, 1, N-1, N, N+1];
    k_nghb1 = k_nghb1(k_nghb1 > 0 & k_nghb1 < N*M+1);
    num_nghb_MA(k) = nnz(ismember(k_nghb1,id_MA2));
end

%--------------------------------------------------------------------------
% Seperate branches using the number of connections
%--------------------------------------------------------------------------
nBranch = 0;
ID = [];
Conn = [];
XY = [];
for i = 1 : length(id_MA2)
    k = id_MA2(i);
%         k = N*(j-1) + i;
        k_nghb1 = k + [-N-1, -N, -N+1, -1, 1, N-1, N, N+1];
        k_nghb2 = k_nghb1 + [-N-1, -N, -N+1, -1, 1, N-1, N, N+1]';
        
        %------------------------------------------------------------------
        % Skip if the point is neither a tip nor a joint of branch
        %------------------------------------------------------------------
        if id_joint(k) == 0 && num_nghb_MA(k) ~= 1
            continue;
        end
        
        %------------------------------------------------------------------
        % Check level-1 neighbor grid points
        %------------------------------------------------------------------
        k_nghb1 = k + [-N-1, -N, -N+1, -1, 1, N-1, N, N+1];
        is_MA = ismember(k_nghb1,id_MA2);
%         is_MA = id_MA_thinned(k_nghb1);
        k_nghb1 = k_nghb1(is_MA);
        
        %------------------------------------------------------------------
        % Find joint to joint connection
        %------------------------------------------------------------------
        id = find(id_joint(k_nghb1));
        if id_joint(k) == 1
            k_nghb1(id) = [];
        elseif any(id)
            for ii = 1 : length(id)
                nBranch = nBranch + 1;
                ID{nBranch} = [k;k_nghb1(id(ii))];
                Conn(nBranch,:) = sort([k,k_nghb1(id(ii))],2);
            end
            k_nghb1(id) = [];
        else
            0;
        end
        
        %------------------------------------------------------------------
        % Loop over the number of connected points
        %------------------------------------------------------------------
        for ik = 1 : length(k_nghb1)
            %--------------------------------------------------------------
            % Collect the connected MA points
            %--------------------------------------------------------------
            k_connection_old = k_nghb1(ik);
            while 1
                k_nghb2 = k_connection_old + [-N-1, -N, -N+1, -1, 1, N-1, N, N+1];
                
                % Remove index out of domain
                [r,c] = ind2sub([N,M],k_nghb2);
                id_boundary = r < 1 | r > N | c < 1 | c > M;
                k_nghb2(id_boundary) = [];
                
                % Remove non-MA neighbor points
                is_MA = id_MA(k_nghb2);
                k_nghb2 = k_nghb2(is_MA);
                k_nghb2 = setdiff(k_nghb2,[k,k_nghb1]);
                
                % Follow joint if multiple neighbor points found
                if any(id_joint(k_nghb2))
                    k_nghb2 = k_nghb2(id_joint(k_nghb2));
                end
                
                % Update connection
                k_connection = vertcat(k_connection_old(:),k_nghb2(:));
                k_connection = unique(k_connection,'stable');
                
                % Quit the loop if there is no update
                if isequal(k_connection_old,k_connection)
                    break;
                end
                
                if any(id_joint(k_connection) == 1) || any(num_nghb_MA(k_connection) == 1)
                    break;
                end
                
                k_connection_old = k_connection;
            end % collecting the MA points connection
            
            %--------------------------------------------------------------
            % Add start point
            %--------------------------------------------------------------
            k_connection = [k;k_connection];
            
            %--------------------------------------------------------------
            % Store the branch ids
            %--------------------------------------------------------------
            nBranch = nBranch + 1;
            ID{nBranch} = k_connection;
            Conn(nBranch,:) = [k_connection(1),k_connection(end)];
        end
end

%--------------------------------------------------------------------------
% Remove redundant connectivity list
%--------------------------------------------------------------------------
Conn = sort(Conn,2);
[~,id] = unique(Conn,'rows');
ID = ID(id);
Conn = Conn(id,:);

% %--------------------------------------------------------------------------
% % Store XY coordinates
% %--------------------------------------------------------------------------
% for i = 1 : length(id)
%     XY{i} = [BX(ID{i}),BY(ID{i})];
% end

%--------------------------------------------------------------------------
% Store the data
%--------------------------------------------------------------------------
MA.nBranch = length(ID);
MA.Size = [N,M];
MA.ID = ID;
MA.ID_Joint = id_joint;
MA.Conn = Conn;
% MA.XY = XY;





