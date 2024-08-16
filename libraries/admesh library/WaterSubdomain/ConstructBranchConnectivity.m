function BranchConnectivity = ConstructBranchConnectivity(JointConnectivity)

%--------------------------------------------------------------------------
% Construct branch connectivity list
%--------------------------------------------------------------------------
BranchConnectivity = cell(size(JointConnectivity,1),1);
for i = 1 : size(JointConnectivity,1)
    BranchConnectivity{i} = [];
    J1 = JointConnectivity(i,1);
    J2 = JointConnectivity(i,2);
    
    id = find(JointConnectivity(:,1) == J1 | JointConnectivity(:,2) == J1);
    id = setdiff(id,i);
    BranchConnectivity{i} = [BranchConnectivity{i}; -id];
    
    id = find(JointConnectivity(:,1) == J2 | JointConnectivity(:,2) == J2);
    id = setdiff(id,i);
    BranchConnectivity{i} = [BranchConnectivity{i}; id];
end