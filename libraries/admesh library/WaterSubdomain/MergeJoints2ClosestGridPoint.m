function MA_new = MergeJoints2ClosestGridPoint(MA)

N = MA.Size(1);
M = MA.Size(2);
ID_Joint = MA.ID_Joint;
ID = MA.ID;
Conn = MA.Conn;
%--------------------------------------------------------------------------
% Find and store redunatant joints data
%--------------------------------------------------------------------------
nDupJoints = 0;
DupJoints = [];
for i = 1 : length(ID)
    for j = 1 : 2
        k = ID{i}(j);
        k_nghb = k + [-N-1, -N, -N+1, -1, 1, N-1, N, N+1];
        
        % Check if it's a joint and has other adjacent joints
        if ID_Joint(k) == 0 || nnz(ID_Joint(k_nghb)) == 0
            continue;
        end
        
        % Collect adjacent joints
        k_joint_connection_old = k;
        while 1
            kk_nghb = k_joint_connection_old + [-N-1, -N, -N+1, -1, 1, N-1, N, N+1];
            id_joint = ID_Joint(kk_nghb);
            k_nghb_joint = kk_nghb(id_joint);
            
            k_joint_connection = vertcat(k_joint_connection_old(:),k_nghb_joint(:));
            k_joint_connection = unique(k_joint_connection);
            if isequal(k_joint_connection_old,k_joint_connection)
                break;
            end
            k_joint_connection_old = k_joint_connection;
        end
        
        % Remove duplicated joints if it's already connected to some branch
%         for ii = 1 : length(k_joint_connection)
%             if nnz(ismember(Conn,k_joint_connection(ii))) > 1
%                 k_joint_connection(ii) = 0;
%             end
%         end
%         k_joint_connection = nonzeros(k_joint_connection);
        
        if isempty(k_joint_connection)
            continue;
        end
        [I,J] = ind2sub([N,M],k_joint_connection);
        mI = mean(I);
        mJ = mean(J);
        dist = sqrt((I - mI).^2 + (J - mJ).^2);
        
        [~,id_closest] = min(dist);
        id_closest = k_joint_connection(id_closest);
        
%         id_closest = sub2ind([N,M],round(mI),round(mJ)); % choose nearest grid point even if it is NOT a joint
        
        nDupJoints = nDupJoints + 1;
        DupJoints.ID_conn{nDupJoints} = k_joint_connection;
        DupJoints.ID_closest(nDupJoints) = id_closest;
    end
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
MA_new = MA;
for i = 1 : nDupJoints
    id_closest = DupJoints.ID_closest(i);
    [I,J] = ind2sub([N,M],id_closest);
    for j = 1 : length(DupJoints.ID_conn{i})
        iconn = DupJoints.ID_conn{i}(j);
        if iconn == id_closest
            continue;
        end
        
        [iBranchs,iEnds] = find(Conn == iconn);
        
        for k = 1 : length(iBranchs)
            iBranch = iBranchs(k);
            iEnd = iEnds(k);
            [II,JJ] = ind2sub([N,M],iconn);
            
            num_add = 0;
            I_add = [];
            J_add = [];
            while 1
                if II == I && JJ == J
                    break;
                end
                
                num_add = num_add+1;
%                 I_add(nadd) = II + sign(I - II);
%                 J_add(nadd) = JJ + sign(J - JJ);
                
                if II < I
                    I_add(num_add) = II + 1;
                elseif II > I
                    I_add(num_add) = II - 1;
                else
                    I_add(num_add) = II;
                end
                
                if JJ < J
                    J_add(num_add) = JJ + 1;
                elseif JJ > J
                    J_add(num_add) = JJ - 1;
                else
                    J_add(num_add) = JJ;
                end
                

                II = I_add(num_add);
                JJ = J_add(num_add);
            end
            

            k_add = sub2ind([N,M],I_add,J_add);
            k_add = k_add(:);
            if abs(I - 1748) < 3 && abs(J - 1762) < 3
                0;
            end
            if iEnd == 1
                k_add = flipud(k_add);
                MA_new.ID{iBranch} = vertcat(k_add,MA_new.ID{iBranch});
            else
                MA_new.ID{iBranch} = vertcat(MA_new.ID{iBranch},k_add);
            end
            MA_new.Conn(iBranch,iEnd) = id_closest;
        end
    end
end

for i = 1 : length(MA_new.ID)
    MA_new.ID{i} = unique(MA_new.ID{i},'stable');
end


% Check duplicated branches





