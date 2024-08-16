function MA = ConstructMedialAxisBranch(id_MA,UIFigure)
%==========================================================================
% ConstructMedialAxisBranch_v2
% - This function does this, this, and this...
% 
% Update history
% ????-??-?? (v1) Written by Younghun Kang
% 2021-03-15 (v2) Younghun: Use built-in 'bwmorph' with 'branchpoints' option to find branch points
% 
%==========================================================================

msg = 'Constructing medial axis branches...';
progdlg = uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');

id_MA = full(id_MA);
N = size(id_MA,1);
M = size(id_MA,2);

%--------------------------------------------------------------------------
% Identify joints based on the number of connected MA points
%--------------------------------------------------------------------------
% JointID = bwmorph(id_MA,'branchpoints');

%--------------------------------------------------------------------------
% Identify joints based on the number of connected MA points
% (Copied & modified a part of "algbwmorph" function for efficient memory
% usage
%--------------------------------------------------------------------------
am_C = bwlookup(id_MA, images.internal.lutbranchpoints());
am_B = bwlookup(id_MA, images.internal.lutbackcount4());
am_C = sparse(am_C);
am_B = sparse(am_B);
am_E = (am_B == 1);
am_FC = ~am_E & am_C; 
clear am_C;

am_Vp = ((am_B == 2) & ~am_E);
am_Vq = ((am_B > 2) & ~am_E);
num_nghb_MA = am_B;
clear am_B am_E;

am_D = bwlookup(full(am_Vq), images.internal.lutdilate());
am_M = (am_FC & am_Vp) & am_D;
clear am_D am_Vp am_Vq;

JointID =  am_FC & ~am_M;
clear am_M am_FC;

%--------------------------------------------------------------------------
% Compute the number of connections
%--------------------------------------------------------------------------
id_MA2 = find(id_MA);
% num_nghb_MA = zeros(length(id_MA2),1);
% for i = 1 : length(id_MA2)
%     k = id_MA2(i);
%     k_nghb1 = k + [-N-1, -N, -N+1, -1, 1, N-1, N, N+1];
%     k_nghb1 = k_nghb1(k_nghb1 > 0 & k_nghb1 < N*M+1);
%     num_nghb_MA(k) = nnz(ismember(k_nghb1,id_MA2));
% end

%--------------------------------------------------------------------------
% Seperate branches using the number of connections
%--------------------------------------------------------------------------
nBranch = 0;
BranchNodes = [];
JointConnectivity = [];
XY = [];
for i = 1 : length(id_MA2)
    k = id_MA2(i);
    
    %------------------------------------------------------------------
    % Skip if the point is neither a tip nor a joint of branch
    %------------------------------------------------------------------
    if JointID(k) == 0 && num_nghb_MA(k) ~= 1
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
    % Find joint-to-joint/tip-to-joint connection
    %------------------------------------------------------------------
    id_k_nghb1_joint = find(JointID(k_nghb1));
    if num_nghb_MA(k) == 1 && ~isempty(id_k_nghb1_joint)
        if length(id_k_nghb1_joint) > 1
            warning(['Wrong detection of tip points are found. ',...
                'It normally happens for short branch (with 1-pixel length), '...
                'So it is removed and proceeds to next branch']);
            continue;
        end
            nBranch = nBranch + 1;
            BranchNodes{nBranch} = [k;k_nghb1(id_k_nghb1_joint)];
            JointConnectivity(nBranch,:) = sort([k,k_nghb1(id_k_nghb1_joint)],2);
        continue;
    elseif JointID(k) == 1 && ~isempty(id_k_nghb1_joint)
        for ii = 1 : length(id_k_nghb1_joint)
            nBranch = nBranch + 1;
            BranchNodes{nBranch} = [k;k_nghb1(id_k_nghb1_joint(ii))];
            JointConnectivity(nBranch,:) = sort([k,k_nghb1(id_k_nghb1_joint(ii))],2);
        end
        k_nghb1(id_k_nghb1_joint) = [];
%         continue;
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
            % (follow closest joint if multiple joints are found)
            k_nghb2_joint = k_nghb2(JointID(k_nghb2));
            if ~isempty(k_nghb2_joint)
                [I1,J1] = ind2sub([N,M],k_connection_old);
                %----------------------------------------------------------
                % 2021-04-12 Younghun:
                % Convert "k_nghb2_joint" to vertical vector in case multiple joints... 
                % !!! Need to check if it should not happen or okay !!!
                %----------------------------------------------------------
                [I2,J2] = ind2sub([N,M],k_nghb2_joint(:)); 
                [~,dist] = knnsearch([I1,J1],[I2,J2]);
                [~,id_k_nghb2_joint] = min(dist);
                k_nghb2 = k_nghb2_joint(id_k_nghb2_joint);
            end
            
            % Update connection
            k_connection = vertcat(k_connection_old(:),k_nghb2(:));
            k_connection = unique(k_connection,'stable');
            
            % Quit the loop if there is no update
            if isequal(k_connection_old,k_connection)
                break;
            end
            
            if any(JointID(k_connection) == 1) || any(num_nghb_MA(k_connection) == 1)
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
        BranchNodes{nBranch} = k_connection;
        JointConnectivity(nBranch,:) = [k_connection(1),k_connection(end)];
    end
    progdlg.Indeterminate = 'off';
    progdlg.Value = i/length(id_MA2);
end

%--------------------------------------------------------------------------
% Remove redundant connectivity list
%--------------------------------------------------------------------------
JointConnectivity = sort(JointConnectivity,2);
[~,id] = unique(JointConnectivity,'rows');
BranchNodes = BranchNodes(id);
JointConnectivity = JointConnectivity(id,:);

% %--------------------------------------------------------------------------
% % Store XY coordinates
% %--------------------------------------------------------------------------
% for i = 1 : length(id)
%     XY{i} = [BX(ID{i}),BY(ID{i})];
% end


%--------------------------------------------------------------------------
% Store the data
%--------------------------------------------------------------------------
if isempty(BranchNodes)
    BranchNodes = {};
end
MA.Size = [N,M];
MA.BranchNodes = BranchNodes;
% MA.XY = XY;





