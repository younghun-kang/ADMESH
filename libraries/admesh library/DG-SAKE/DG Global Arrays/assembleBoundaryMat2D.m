function [B,IL,IR,BC] = assembleBoundaryMat2D(MESH,b,Ae,L,ndof,nbp,nelems)
% *************************************************************************
% Assembly/Compute the global advection matrices
%
% Modified from Ethan Kubatko's program
%
% Author:
% Dustin West
% Assistant Engineer | Hazen and Sawyer               
% dwest@hazenandsawyer.com | hazenandsawyer.com
% *************************************************************************

%------------------------------------------------------------------
% Initialize variables
%------------------------------------------------------------------
E           = MESH.edges;               % edge list
EA          = MESH.edgeAttachments;     % edge attachments
ET          = MESH.edgeType;            % edge type

nedges      = size(E,1);                % # of edges
nfb         = sum(ET == 0 | ET == 1);   % # of free boundary edges
nce         = sum(ET == -2);            % # of converging edges
%nce         = sum(ET == -2 | ET == -3 );% # of converging edges

% Initialize global boundary flux matrix and sparse matrix indexing
[B.ii,ib,jb] = deal(zeros(ndof*(nedges - nfb/2)*2*nbp ,1));

% Initialize "left", "right" indexing & BC vector
[IL.ii,IR.ii,BC.ii]  = deal(zeros((nedges + nce)*nbp,1));

% Initialize vector indexing for IL & IR
lv = 1:nbp;

% Initialize column subscript indexing for global matrices
ic = ones(ndof,1)*(1:nbp);

% Initialize indexing vector for i & j
nId = ndof*nbp; id = 1:nId;

%------------------------------------------------------------------
% Loop over each edge
%------------------------------------------------------------------
for j = 1:nedges
    
    %--------------------------------------------------------------
    % LEFT ELEMENT ASSEMBLY
    %--------------------------------------------------------------
    
    % Retrieve "left" element attached to edge
    iL = EA{j}(1);
    
    % Determine local edge number of edge j from within element iL
    kL = getLocalEdgeNumber(MESH,iL,E(j,:));
    
    % Place the element boundary matrix in the correct global slot
    B.ii(id) = L{kL}(iL)/Ae(iL)*b.ii{kL};
    
    % Populate sparse matrix indexing
    ib(id) = (1:ndof)'*ones(1,nbp) + ndof*(iL-1);
    jb(id) = ic;
    
    %--------------------------------------------------------------
    % RIGHT ELEMENT ASSEMBLY
    %--------------------------------------------------------------
    
    if length(EA{j}) == 2 % Perform "right" element assembly
        
        % Retrieve "right" element attached to edge
        iR = EA{j}(2);
        
        % Determine local edge number of edge j from within element iR
        kR = getLocalEdgeNumber(MESH,iR,E(j,:));
        
        % Set up "right" element assembly based on edge type
        switch ET(j)
            
            case -1 % FLOW EDGE
                
                % Increment index for boundary matrix
                id = id + nId;
                
                % Place the element boundary matrix in the correct global slot
                B.ii(id) = -L{kR}(iR)/Ae(iL)*fliplr(b.ii{kR});
                
                % Populate sparse matrix indexing
                ib(id) = (1:ndof)'*ones(1,nbp) + ndof*(iR-1);
                jb(id) = ic;
                
                % Determine the edge to global integration point index
                aL = ((iL-1)*3*nbp + (kL-1)*nbp + 1);
                aR = ((iR-1)*3*nbp + (kR-1)*nbp + 1);
                
                % Store "left" & "right" index
                IL.ii(lv)   = aL:aL+nbp-1;
                IR.ii(lv)   = fliplr(aR:aR+nbp-1);
                
                % Store internal neighbor boundary condition
                BC.ii(lv) = -1;
                
                % Increase index for IL, IR and BC
                lv = lv + nbp;
                
            case -3 % FLOW EDGE (1D routing edge)
                
                % Increment index for boundary matrix
                id = id + nId;
                
                % Place the element boundary matrix in the correct global slot
                B.ii(id) = -L{kR}(iR)/Ae(iL)*fliplr(b.ii{kR});
                
                % Populate sparse matrix indexing
                ib(id) = (1:ndof)'*ones(1,nbp) + ndof*(iR-1);
                jb(id) = ic;
                
                % Determine the edge to global integration point index
                aL = ((iL-1)*3*nbp + (kL-1)*nbp + 1);
                aR = ((iR-1)*3*nbp + (kR-1)*nbp + 1);
                
                % Store "left" & "right" index
                IL.ii(lv)   = aL:aL+nbp-1;
                IR.ii(lv)   = fliplr(aR:aR+nbp-1);
                
                % Store internal neighbor boundary condition
                BC.ii(lv) = -1;
                
                % Increase index for IL, IR and BC
                lv = lv + nbp;
                
            case 0 % NO FLOW EDGE
                
                % Determine the edge to global integration point index
                aL = ((iL-1)*3*nbp + (kL-1)*nbp + 1);
                
                % Store "left" & "right" index & edge id
                IL.ii(lv)   = aL:aL+nbp-1;
                IR.ii(lv)   = IL.ii(lv);
                
                % Copy edge boundary type
                BC.ii(lv) = 0;
                
                % Increase index for IL, IR and BC
                lv = lv + nbp;
                
            case -2 % CASCADE EDGE
                
                % complete "left" element
                
                % Determine the 2D edge to global integration point index
                aL = ((iL-1)*3*nbp + (kL-1)*nbp + 1);
                
                % Store 2D "left" & "right" index & edge id
                IL.ii(lv)   = aL:aL+nbp-1;
                IR.ii(lv)   = IL.ii(lv);
                
                % Copy edge boundary type
                BC.ii(lv) = -2;
                
                % Increase index for IL & IR
                lv = lv + nbp;
                
                % Repeat for "right" element and treat as a "left" element
                
                % Increment id
                id = id + nId;
                
                % Increment column index for next edge
                ic = ic + nbp;
                
                % Place the element boundary matrix in the correct global slot
                B.ii(id) = L{kR}(iR)/Ae(iR)*b.ii{kR};
                
                % Populate sparse matrix indexing
                ib(id) = (1:ndof)'*ones(1,nbp) + ndof*(iR-1);
                jb(id) = ic;
                
                % Determine the edge to global integration point index
                aR = ((iR-1)*3*nbp + (kR-1)*nbp + 1);
                
                % Store "left" & "right" index & edge id
                IL.ii(lv)   = aR:aR+nbp-1;
                IR.ii(lv)   = IL.ii(lv);
                
                % Store internal neighbor boundary condition
                BC.ii(lv) = -2;
                
                % Increase index for IL & IR
                lv = lv + nbp;
                
        end
        
    else % FREE BOUNDARY EDGE
        
        % Determine the edge to global integration point index
        aL = ((iL-1)*3*nbp + (kL-1)*nbp + 1);
        
        % Store "left" & "right" index & edge id
        IL.ii(lv)   = aL:aL+nbp-1;
        IR.ii(lv)   = IL.ii(lv);
        
        % Copy edge boundary type
        BC.ii(lv) = ET(j);
        
        % Increase index for IL & IR
        lv = lv + nbp;
        
    end
    
    % Increment id
    id = id + nId;
    
    % Increment column index for next edge
    ic = ic + nbp;
    
end

% Construct sparse boundary matrix
B.ii   = sparse(ib,jb,B.ii,nelems*ndof,(nedges + nce)*nbp);

%=================================================================
% % SUBROUTINES
%=================================================================
    function k = getLocalEdgeNumber(MESH,ix,E)
        
        k1 = find(MESH.ConnectivityList(ix,:)==E(1,1));
        k2 = find(MESH.ConnectivityList(ix,:)==E(1,2));
        k  = abs(k1-k2) + mod(max(k1,k2),3);
        
    end

end
