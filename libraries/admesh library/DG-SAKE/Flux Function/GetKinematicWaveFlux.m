function [F,G,ck,S] = GetKinematicWaveFlux(N,X,Y,ZL,ZR,MESH,n)
%************************************************************************
%************************************************************************
% Compute TIN Kinematic Wave Properties
%
% Input: MESH: Triangular Mesh
%        Fz:   Elevation Function
%        niElem: Number of element integration points in all elements
%        niEdge: Number of eedge integration points in all elements
%
% Output: Sx: Global slope in x direction
%         Sy: Global slope in y direction
%         So: Global slope magnitude
%          f: Flux function in x
%          g: Flux function in y
%
% Author:
% Dustin West
% Assistant Engineer | Hazen and Sawyer               
% dwest@hazenandsawyer.com | hazenandsawyer.com
%************************************************************************
%************************************************************************

%========================================================================
% Assemble 2D Flux Components
%========================================================================

m  = 5/3;   % Turbulent Flow Conditions
%K  = 1.49;  % When using english units
K = 1   ;  % When using SI units
   
%------------------------------------------------------------------------
% Assemble Global Slope Components
%------------------------------------------------------------------------
nelems      = size(MESH.ConnectivityList,1);    % # of elements
nap         = numel(X.elem2)/nelems;            % # of element points
nbp         = numel(X.edge2)/nelems;            % # of edge points

% Compute slope in x and in y
[~,dzdx,dzdy] = triSlopeVector(MESH); mz = sqrt(dzdx.^2 + dzdy.^2);

% Initialize Slope Component vectors for element points
[Sx.elem2,Sy.elem2,So.elem2] = deal(zeros(nap*nelems,1));

% Copy element slope to all element integration points
Sx.elem2(:) = ones(nap,1)*dzdx'; % Slope in x at points
Sy.elem2(:) = ones(nap,1)*dzdy'; % Slope in y at points
So.elem2(:) = ones(nap,1)*mz';   % Magnitude


% Initialize Slope Component vectors for edge points
[Sx.edge2,Sy.edge2,So.edge2] = deal(zeros(nbp*nelems,1));

Sx.edge2(:) = ones(nbp,1)*dzdx'; % Slope in x at points
Sy.edge2(:) = ones(nbp,1)*dzdy'; % Slope in y at points
So.edge2(:) = ones(nbp,1)*mz';   % Magnitude

%------------------------------------------------------------------------
% Create global flux function for element integration points
%------------------------------------------------------------------------

% Flux function for flow rate in x
alpha   = (K./n.elem2).*(Sx.elem2./sqrt(So.elem2)); % element
F.elem2 = @(h)( alpha.*(h.^m) );         % element

% Flux function for flow rate in y
alpha   = (K./n.elem2).*(Sy.elem2./sqrt(So.elem2)); % element
G.elem2 = @(h)( alpha.*(h.^m) );         % element

alpha   = (K./n.edge2).*(Sx.edge2./sqrt(So.edge2)); % edge
F.edge2 = @(h)( alpha.*(h.^m) );         % edge

alpha   = (K./n.edge2).*(Sy.edge2./sqrt(So.edge2)); % edge
G.edge2 = @(h)( alpha.*(h.^m) );         % edge

%------------------------------------------------------------------------
% Create kinemeatic wave celerity function for propogation speed
%------------------------------------------------------------------------
% Derivative of Flux function for flow rate in x
alpha   = (K.*m./n.edge2).*(Sx.edge2./sqrt(So.edge2)); % edge
cx      = @(h,id)( alpha(id).*(h(id).^(m-1)) );

alpha   = (K.*m./n.edge2).*(Sy.edge2./sqrt(So.edge2)); % edge
cy      = @(h,id)( alpha(id).*(h(id).^(m-1)) );

% Compute anonymous propogation speed function
ck.ii  = @(h,int,ext)...
    max(...
    abs(cx(h,int).*N.x(int) + cy(h,int).*N.y(int) ),...
    abs(cx(h,ext).*N.x(int) + cy(h,ext).*N.y(int) ));

%========================================================================
% Assemble 1D Flux Components
%========================================================================

%------------------------------------------------------------------
% Assemble 1D?
%------------------------------------------------------------------
if ~isfield(MESH, 'channelNetwork'); S = nan; return; end

%------------------------------------------------------------------
% Compute slope components for each element
%------------------------------------------------------------------
nelems  = sum(cellfun(@(x) size(x,1),MESH.channelNetwork)-1);   % # of elements
ngauss  = size(X.elem1,1);                               % Number of gauss points

%------------------------------------------------------------------
% Compute element slopes and assign channel widths
%------------------------------------------------------------------
                     
% Initialize vectors
S.i         = zeros(ngauss*nelems,1);
F.elem1.it  = true(ngauss,nelems);
F.elem1.ir  = false(ngauss,nelems);

% Initialize indexing vector
ix = 0;

% Loop over each cascade segment in MESH.cascadeNetwork
for k = 1:length(MESH.channelNetwork)
    
    % Assign nodes
    nodes = MESH.channelNetwork{k}(:,1);
    
    % Get (x,y) coordinates
    x = MESH.Points(nodes,1);
    y = MESH.Points(nodes,2);
    z = MESH.Points(nodes,3);
    
    % Compute element sizes in string
    h = sqrt(diff(x).^2 + diff(y).^2);
    
    % Compute difference in elevation
    dz = diff(z);
    
    % compute indexing vector
    ix = (ix(end)+1):(ix(end)+length(dz)*ngauss);
    
    % Compute slope of each element
    S.i(ix) = ones(ngauss,1)*abs(dz./h)';
    
end

% Check for near zero slopes
%tol = 1*10^-(4.5);
tol = eps;
if any(S.i <= tol)
    disp(' There is a 1D element with less then or zero slope')
    S.i(S.i <= tol) = mean(S.i);
end

n.elem1 = reshape(n.elem1,size(X.elem1,1),size(X.elem1,2));

% Designate element cross section type
for k = 1:length(MESH.Constraints)
    
    if MESH.Constraints(k).num == 18
        
        % Get node string
        nodes = MESH.Constraints(k).nodes;
        
        % Get element edges
         x = [MESH.Points(nodes(1:end-1),1),MESH.Points(nodes(2:end),1)]' ;
         %y = [MESH.Points(nodes(1:end-1),2),MESH.Points(nodes(2:end),2)]';
        
         for j = 1:size(x,2)
             
             elemi = find((x(1,j) == X.elem1(1,:) & x(2,j) == X.elem1(end,:)) | (x(2,j) == X.elem1(1,:) & x(1,j) == X.elem1(end,:)) );
             
             if length(elemi) > 1
                 error('Hey!')
             end

             % Store indices
             F.elem1.ir(:,elemi) = true;
             F.elem1.it(:,elemi) = false;
             
             % Set mannings n to 0.03
             n.elem1(:,elemi) = 0.035;
             
         end
        
    end
    
end

n.elem1    = n.elem1(:);
F.elem1.ir = F.elem1.ir(:);
F.elem1.it = F.elem1.it(:);

end
