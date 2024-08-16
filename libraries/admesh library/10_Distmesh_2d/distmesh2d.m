function MESH = distmesh2d(PTS,phi,MeshFun,xyzFun,hmin,Settings,UIFigure,pH,AdvSettings)
% distmesh2d - Generates a mesh based on mesh size h
%
% Syntax:  [p,t] = distmesh2d(DistanceFun,MeshSizeFun,hmin,guiFig)
%
% Inputs:
%    DistanceFun - gridded interpolant of the Distance Function
%    MeshSizeFun - gridded interpolant of the Mesh Size Function
%    hmin        - minimum element size
%    guiFig - handle that identifies the figure
%
% Outputs:
%    p - points of delaunay triangulation
%    t - connectivity list
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

% DISTMESH2D 2-D Mesh Generator using Distance Functions.
% distmesh2d.m v1.1
% Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.
%
% Author of Adjustments: Dustin West, Colton Conroy, Ethan Kubatko
% The Ohio State University
% email address: dww.425@gmail.com
% August 2013; Last revision: 08-August-2013
%--------------------------- BEGIN CODE -----------------------------------

%--------------------------------------------------------------------------
% GUI check - Check if user wants to view mesh generation
%--------------------------------------------------------------------------
viewMesh = strcmpi(Settings.View.Status,'On'); axes(pH);

%--------------------------------------------------------------------------
% Initialize variables
%--------------------------------------------------------------------------
% Specifies how far the points can move (relatively) before a retriangulation 
ttol = 0.1;
% "Internal pressure" of truss
Fscale  = 1.2;    
% Time step in Euler's method
deltat  = .2;     
% Tolerance in the geometry evaluations.
geps    = .01*hmin; 
% Frequency in which to check for nodes that are too close to one another
densityctrlfreq= 50;  
% Number of iterations
niter   = 5*1100;  
niter  = 100;
% Initialize current positions
pold=inf; 
% Intitalize current element quality
qold    = 0; 

if ~isempty(AdvSettings)
    ttol       = AdvSettings.ttol;
    Fscale     = AdvSettings.Fscale;
    deltat     = AdvSettings.deltat;
    geps_scale = AdvSettings.geps_scale;
    niter      = AdvSettings.niter;
    qold       = AdvSettings.qold;

    geps = geps_scale*hmin;
end

msg = 'Creating initial distribution of nodal points...';
progdlg = uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');

%--------------------------------------------------------------------------
% Create initial distribution in bounding box (equilateral triangles)
%--------------------------------------------------------------------------
p = createInitialPointList(PTS,hmin);

%--------------------------------------------------------------------------
% Remove points outside the region, apply the rejection method
%--------------------------------------------------------------------------
p = rejectionMethod(p,phi.f,MeshFun,geps);

%**************************************************************************
% Younghun temporarly added
%**************************************************************************

ConstraintsXY = [];
if ~isempty(PTS.Constraints)
for i = 1 : length(PTS.Constraints)
    ConstraintsXY = vertcat(ConstraintsXY, PTS.Constraints(i).xy, [nan nan]);
end
ConstraintsXY1 = vertcat(PTS.Constraints(:).xy);
end
% PTS.Constraints = []; 

%--------------------------------------------------------------------------
% Apply mesh constraints and concatenate with p if constraints exist.
%--------------------------------------------------------------------------
[p,nC,C,MESH] = GetMeshConstraints(p,hmin,PTS);
N = size(p,1); % number of nodes

in = ((nC+1):N)'; % Vector of non-pfix indices

msg = 'Generating mesh...';
progdlg = uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg);
progdlg.Value = 1/niter;

drawnow;

%--------------------------------------------------------------------------
% Remove previous plot before starting
%--------------------------------------------------------------------------

if viewMesh == 1; 
    
    h = findobj(pH,'tag', 'Edge Structure');
    h = [h;findobj(pH,'tag', 'Internal Constraint')];
    h = [h;findobj(pH,'tag', 'External Constraint')];
    h = [h;findobj(pH,'tag', 'Open Ocean')];
    h = [h;findobj(pH,'tag', 'Line Constraint')];
    
    if ~isempty(h); delete(h); end
    
    h = findobj('tag', 'Mesh');
    
    if ~isempty(h); delete(h); end 

end

Ftot = inf;
%--------------------------------------------------------------------------
% Generate Mesh
%--------------------------------------------------------------------------
for k = 1:niter
    
    drawnow; % for graphics

    % Retriangulation by the Delaunay algorithm 
    if ( max(deltat*sqrt(sum(Ftot.^2,2)))/hmin > ttol )  % Any large movement?

        % Save current positions
        pold=p;
        
        % Perform Delaunay Triangulation
        if isempty(C)
            dt = delaunayTriangulation(p);
        else
            dt = delaunayTriangulation(p,C);
        end
        
        %------------------------------------------------------------------
        % 2021-03-18 Younghun: fix a issue (not sure when it happens) by a temporary solution
        %------------------------------------------------------------------
        if ~isequal(p,dt.Points)
            warning([...
                'The routine ''delaunayTriangulation'' results more/less points.',...
                'As a temporal solution, modify the variable dt so that match to original p.']);
            %                 id_p = (~ismember(dt.Points,p,'rows'));
            %                 dt.Points = p;
            
            %                 warning([...
            %                     'The routine ''delaunayTriangulation'' results more/less points.',...
            %                     'As a temporal solution, replace original p with Points from ''delaunayTriangulation''.']);
            p = dt.Points;
            pold = p;
            t = dt.ConnectivityList;
            N = size(p,1);
            in = ((nC+1):N)'; % Vector of non-pfix indices
            % id = find(ismember(dt.ConnectivityList(:,1),find(id_p)));
            % I = ind2sub(size(dt.ConnectivityList),id);
            % dt.ConnectivityList(I,:) = [];
            
        end

        % Compute centroids & Interpolate distances 
        ind = phi.f((p(dt(:,1),:)+p(dt(:,2),:)+p(dt(:,3),:))/3) < -geps;
        
        % Keep interior triangles
        t = sort(dt(ind,:),2);
        

        % Describe each bar by a unique pair of nodes
        bars = unique([t(:,[1,2]);t(:,[1,3]);t(:,[2,3])],'rows');

        % Graphical output of the current mesh
        if viewMesh
            
            % Delete current plot
            h = findobj('tag', 'Mesh'); if ~isempty(h); delete(h); end
            cla(pH);
            % Plot mesh
            patch(pH,...
                'vertices',p,...
                'faces',t,...
                'edgecol',[0 .4 .8],...
                'facecol','none',...
                'Tag','Mesh');
            
            
            if ~isempty(ConstraintsXY)
                plot(pH,ConstraintsXY(:,1),ConstraintsXY(:,2),':','color',[0.5 0 0],'Tag','Mesh');
            end
            
            % Plot Constraints
            if ~isempty(C)
     
                patch(pH,'vertices',p,'faces',C,'edgecol','r','facecol','none','Tag','Mesh');
                
            end

        end
        
    end
        
    % Move mesh points based on bar lengths L and forces F
    barvec  = p(bars(:,1),:)-p(bars(:,2),:);              % List of bar vectors
    L       = sqrt(sum(barvec.^2,2));                     % L = Bar lengths
    hbars   = MeshFun((p(bars(:,1),:)+p(bars(:,2),:))/2); % hbar = mesh size
    L0      = hbars*Fscale*sqrt(sum(L.^2)/sum(hbars.^2)); % L0 = Desired lengths
    
%     % Density control
%     if( (mod(k,densityctrlfreq) == 0) && (k < niter-5) )
%         
%         set(sb.ProgressBar,'Value',k)
%         sb.setText(['Generating mesh... ' num2str((k/niter)*100,'%.0f') '%'])
%         
%         if any(L0>2*L) % Remove points that are too close together
%             
%             p(setdiff(reshape(bars(L0>2*L,:),[],1),1:nC),:)=[];
%             
%             N=size(p,1); 
%             pold=inf;
%             
%             in = ((nC+1):N)';
%             
%             continue;
%             
%         elseif(k > niter/2)
%             
%             p = BoundaryDensityControl(p,t,C);
%             
%             p = ConstraintDensityControl(p,nC,C,MeshFun);
%             
%             N=size(p,1); 
%             pold=inf; 
%             in = ((nC+1):N)';
%             
%             continue;
%             
%         end
%         
%     end

    F            = max(L0-L,0);         % Bar forces (scalars)
    Fvec         = F./L*[1,1].*barvec;  % Bar forces (x,y components)
    Ftot         = full(sparse(bars(:,[1,1,2,2]),ones(size(F))*[1,2,1,2],[Fvec,-Fvec],N,2));
    Ftot_nC      = Ftot(1:nC,:); % 2021-04-26 Younghun
    Ftot(1:nC,:) = 0;                   % Force = 0 at fixed points
    p            = p+deltat*Ftot;       % Update node positions

    %--------------
    % 2021-04-25 Younghun: Apply projection to the nodes near channels
    %--------------
    AddNodes2Constraints = 1;
    if AddNodes2Constraints && k > niter*0.8 && ~isempty(PTS.Constraints) %&& nAddRemove >= 5
        % Find id of weak constraints
        id_wc = find(-phi.f(p) < hmin*2 & phi.f(p) < 0);
        id_wc = setdiff(id_wc,1:nC);
        cXY = struct2table(PTS.Constraints,'AsArray',true);
        cXY = cXY.xy;

        [Vx,Vy,D] = Compute8SSED_v3(cXY,p(id_wc,1),p(id_wc,2),hmin);
%         [~,~,~,XY1,XY2] = Compute8SSED_v3({p(1:nC,:)},p(id_wc,1),p(id_wc,2),hmin);
        [~,~,~,XY1,XY2] = VectorDistanceTransform(cXY,p(id_wc,:));
        D1 = sqrt(sum((XY1 - XY2).^2,2));
        id_wc_sub = find(abs(D) > 0 & abs(D) < hmin*0.5 &...
            ~isnan(D) & ~isinf(D) & D1 > hmin*2);
        id_wc = id_wc(id_wc_sub);
%         p_add_wc = p(id_wc,:) + [Vx(id_wc_sub), Vy(id_wc_sub)];
        p_add_wc = (XY1+XY2)/2;

        p_add_wc = p_add_wc(id_wc_sub,:);

        
        [~,J1] = (ismember(XY1(id_wc_sub,:),p,'rows'));
        [~,J2] = (ismember(XY2(id_wc_sub,:),p,'rows'));
        id_sc_nearest = [J1,J2];
        
        if nnz(id_sc_nearest == 0) > 0
            error('something is wrong..');
        end
        [~,J1] = ismember(id_sc_nearest(:,[1 2]),C,'rows');
        [~,J2] = ismember(id_sc_nearest(:,[2 1]),C,'rows');
        iRows_C = J1;
        iRows_C(J1 == 0) = J2(J1 == 0);
        I = false(length(id_wc),1);
        I(iRows_C ~= 0) = true;

        % Make only one point can project to one constraint segment at one timestep
        [~,I1] = unique(iRows_C,'stable');
        I(setdiff(1:length(id_wc),I1)) = false;
        
        % Do not project to the current constraint node
        I1 = ~ismember(p_add_wc,p,'rows');
        I = I & I1(:);
        
        d1 = sqrt(sum((p(id_wc,:) - XY1(id_wc_sub,:)).^2,2));
        d2 = sqrt(sum((p(id_wc,:) - XY2(id_wc_sub,:)).^2,2));
        d = min(d1,d2);
        R = abs(D(id_wc_sub))./d;
        
        I_add = I & R < 0.5;
        iRows_C = iRows_C(I_add);
        id_wc = id_wc(I_add);
        id_wc_sub = id_wc_sub(I_add);
        p_add_wc = p_add_wc(I_add,:);
        
        % Add new points to constraint
        for i = 1 : length(iRows_C)
            irC = iRows_C(i);
            C(end+1,:) = [C(irC,1) id_wc(i)];
            C(end+1,:) = [id_wc(i), C(irC,2)];
        end
        C(iRows_C,:) = [];

        p(id_wc,:) = p_add_wc;

        % Reorder node id
        id_freePoints = setdiff(1:size(p,1),1:nC);
        id_freePoints = setdiff(id_freePoints,id_wc);
        id_new = [(1:nC)'; id_wc(:); id_freePoints(:)];
        p = p(id_new,:);
%         pold = pold(id_new,:);

        [~,C] = ismember(C,id_new);

        nC = nC+length(id_wc);
%         pold = p;
        N = size(p,1);
        in = ((nC+1):N)'; % Vector of non-pfix indices

        % Perform Delaunay Triangulation
        if isempty(C)
            dt = delaunayTriangulation(p);
        else
            dt = delaunayTriangulation(p,C);
        end
        p = dt.Points;
        N = size(p,1);
        in = ((nC+1):N)'; % Vector of non-pfix indices
            
        % Compute centroids & Interpolate distances
        ind = phi.f((p(dt(:,1),:)+p(dt(:,2),:)+p(dt(:,3),:))/3) < -geps;
        
        % Keep interior triangles
        t = sort(dt(ind,:),2);

        % Describe each bar by a unique pair of nodes
        bars = unique([t(:,[1,2]);t(:,[1,3]);t(:,[2,3])],'rows');
    end  
    
    %----------------------------------------------------------------------
    % 2021-04-26 Younghun: Project constraints nodes back to channels
    %----------------------------------------------------------------------
    MoveWeakConstraints = 0;
    if MoveWeakConstraints && k > niter*0.7 && ~isempty(PTS.Constraints)
        Ftot_nC2 = sqrt(sum(Ftot_nC.^2,2));
        id = Ftot_nC2 > 0.5*hmin;
        Ftot_nC(id,:) = Ftot_nC(id,:).*(0.5*hmin./Ftot_nC2(id));
        p(1:nC,:) = p(1:nC,:) + deltat*Ftot_nC;
        id_wc = 1:nC;
        cXY = struct2table(PTS.Constraints,'AsArray',true);
        cXY = cXY.xy;
        cXY = vertcat({[PTS.Poly.x(:),PTS.Poly.y(:)]}, cXY(:));
        
        [Vx,Vy,D] = Compute8SSED_v3(cXY,p(id_wc,1),p(id_wc,2),hmin);
        id_wc_sub = abs(D) > 0 & abs(D) < hmin*0.5 & ~isnan(D) & ~isinf(D);
        p(id_wc(id_wc_sub),:) = p(id_wc(id_wc_sub),:) + [Vx(id_wc_sub),Vy(id_wc_sub)];
%         p(id_wc(~id_wc_sub),:) = p(id_wc(~id_wc_sub),:) - deltat*Ftot_nC(~id_wc_sub,:)*0.01;
    end  
    
    % Bring outside points back to the boundary
    p(in,:) = projectBackToBoundary(phi,p(in,:));
    
    %------------------------------------------------------------------
    % 2021-04-13 Younghun: Remove nodes near constraint lines
    %------------------------------------------------------------------
    RemoveNodesNearConstraints = 1;
    if RemoveNodesNearConstraints && k > niter*0.7 && ~isempty(PTS.Constraints)
        id_rm = find(-phi.f(p) < hmin*2);
        id_rm = setdiff(id_rm,1:nC);
        cXY = struct2table(PTS.Constraints,'AsArray',true);
        cXY = cXY.xy;
        % temp = vertcat({[PTS.Poly.x,PTS.Poly.y]}, temp(:));
        
        [~,~,D] = Compute8SSED_v3(cXY,p(id_rm,1),p(id_rm,2),hmin*0.25);
%         [~,~,~,XY1,XY2] = Compute8SSED_v4({p(1:nC,:)},C,p(id_rm,1),p(id_rm,2),hmin);
%         D1 = sqrt(sum((XY1 - XY2).^2,2));
        
        id2 = abs(D) > 0 & abs(D) < hmin*0.5 & ~isnan(D) & ~isinf(D);
        id_rm = id_rm(id2);
        p(id_rm,:) = [];
%         pold(id_rm,:) = [];
        N = size(p,1);
        in = ((nC+1):N)'; % Vector of non-pfix indices
        
        % Perform Delaunay Triangulation
        if isempty(C)
            dt = delaunayTriangulation(p);
        else
            dt = delaunayTriangulation(p,C);
        end
        p = dt.Points;
        N = size(p,1);
        in = ((nC+1):N)'; % Vector of non-pfix indices
        
        % Compute centroids & Interpolate distances
        ind = phi.f((p(dt(:,1),:)+p(dt(:,2),:)+p(dt(:,3),:))/3) < -geps;
        
        % Keep interior triangles
        t = sort(dt(ind,:),2);

        % Describe each bar by a unique pair of nodes
        bars = unique([t(:,[1,2]);t(:,[1,3]);t(:,[2,3])],'rows');
    end

%         p(in,:) = projectBackToBoundary(phi,p(in,:));

%     %----------------------------------------------------------------------
%     % 2021-04-08 Younghun: Apply magnetic-like force
%     %----------------------------------------------------------------------
%     id = find(-phi.f(p) < hmin*2);
%     temp = struct2table(PTS.Constraints);
%     temp = temp.xy;
%     temp = vertcat({[PTS.Poly.x,PTS.Poly.y]}, temp(:));
%     
%     [Vx,Vy,D] = Compute8SSED_v3(temp,p(id,1),p(id,2),hmin);
%     id2 = abs(D) < hmin*0.4 & ~isnan(Vx) & ~isnan(Vy) & ~isinf(Vx) & ~isinf(Vy);
%     p(id(id2),:) = p(id(id2),:) + [Vx(id2),Vy(id2)];
    

   
    %------------------------------------------------------------------
    % 2021-04-15 Younghun: Perform triangulation once again because the updated
    % positions (p) possibly make elements crossing the constraints
    %------------------------------------------------------------------
    if k > niter*0.7
        % Perform Delaunay Triangulation
        if isempty(C)
            dt = delaunayTriangulation(p);
        else
            dt = delaunayTriangulation(p,C);
        end        
%         if ~isequal(p,dt.Points)
%             warning([...
%                 'The routine ''delaunayTriangulation'' results more/less points.',...
%                 'As a temporal solution, modify the variable dt so that match to original p.']);
%             %                 id_p = (~ismember(dt.Points,p,'rows'));
%             %                 dt.Points = p;
%             
%             %                 warning([...
%             %                     'The routine ''delaunayTriangulation'' results more/less points.',...
%             %                     'As a temporal solution, replace original p with Points from ''delaunayTriangulation''.']);
%             p = dt.Points;
%             pold = p;
%             t = dt.ConnectivityList;
%             N = size(p,1);
%             in = ((nC+1):N)';
%         end
        % Compute centroids & Interpolate distances
        if any(max(dt.ConnectivityList) > size(p,1))
            0;
        end
        ind = phi.f((p(dt(:,1),:)+p(dt(:,2),:)+p(dt(:,3),:))/3) < -geps;
        
        % Keep interior triangles
        t = sort(dt(ind,:),2);
    end
    
    % Check element quality. Keep track of best triangulation
    if k > (niter-50)
        [q, ~] = MeshQuality(p,t,0,'Triangle');
        if q > qold; P = p; T = t; qold = q; Csave = C; end
    end
    
    if mod(k/niter,.01) == 0
        progdlg.Value = k/niter;
    end
    
end

close(progdlg);

%--------------------------------------------------------------------------
% Younghun: Re-construct triangulations and remove elements more precisely
%--------------------------------------------------------------------------
if isempty(C)
    T = delaunayTriangulation(P);
else
    T = delaunayTriangulation(P,Csave);
end
P = T.Points;
mP = (P(T(:,1),:)+P(T(:,2),:)+P(T(:,3),:))/3;
ind = PointsInDomain3(mP(:,1),mP(:,2),PTS);
ind = ind & phi.f(mP) < -geps;
% Keep interior triangles
T = sort(T(ind,:),2);
        
% %--------------------------------------------------------------------------
% % 2021-04-13 Younghun: reject nodes near constraint lines
% %--------------------------------------------------------------------------
% id = find(-phi.f(P) < hmin*2);
% temp = struct2table(PTS.Constraints);
% temp = temp.xy;
% % temp = vertcat({[PTS.Poly.x,PTS.Poly.y]}, temp(:));
% 
% [Vx,Vy,D] = Compute8SSED_v3(temp,P(id,1),P(id,2),hmin);
% id2 = abs(D) > 0 & abs(D) < hmin*0.4 & ~isnan(Vx) & ~isnan(Vy) & ~isinf(Vx) & ~isinf(Vy);
% P(id(id2),:) = [];

% %--------------------------------------------------------------------------
% % 2021-03-18 Younghun: Perform triangulation once again because the updated
% % positions (p) possibly make elements crossing the constraints
% %--------------------------------------------------------------------------
% % [~,~,~,MESH] = GetMeshConstraints(p,hmin,PTS);
% if isempty(C)
%     dt = delaunayTriangulation(P);
% else
%     dt = delaunayTriangulation(P,C);
%     % 2021-03-18 Younghun added below to fix a issue (not sure when it happens) by a temporary solution
%     if ~isequal(P,dt.Points)
%         warning([...
%             'The routine ''delaunayTriangulation'' results more/less points.',...
%             'As a temporal solution, modify the variable dt so that match to original p.']);
% %         id_p = (~ismember(dt.Points,p,'rows'));
% % %         dt.Points = P;
%         P = dt.Points;
%         C = dt.Constraints;
%         % id = find(ismember(dt.ConnectivityList(:,1),find(id_p)));
%         % I = ind2sub(size(dt.ConnectivityList),id);
%         % dt.ConnectivityList(I,:) = [];
%         
%     end
% end
% % Compute centroids & Interpolate distances
% ind = phi.f((P(dt(:,1),:)+P(dt(:,2),:)+P(dt(:,3),:))/3) < -geps;
% 
% % Keep interior triangles
% T = sort(dt(ind,:),2);

% %--------------------------------------------------------------------------
% % 2021-04-12 Younghun: Apply magnetic force as a post-process
% %--------------------------------------------------------------------------
% id = find(-phi.f(p) < hmin);
% temp = struct2table(PTS.Constraints);
% temp = temp.xy;
% temp = vertcat({[PTS.Poly.x,PTS.Poly.y]}, temp(:));
% 
% [~,~,D] = Compute8SSED_v3(temp,p(id,1),p(id,2),hmin);
% id2 = D < 0 & abs(D) < hmin/2;
% 
% temp = cell2mat(temp);
% [id2B,dist2B] = knnsearch(temp,p(id,:));
% id2 = abs(D) < hmin/2 & dist2B < hmin;
% % id2 = dist2B < hmin*0.5;
% p1 = p;
% p1(id(id2),:) = temp(id2B(id2),:);


%--------------------------------------------------------------------------
% Clean up 
%--------------------------------------------------------------------------
msg = 'Cleaning up final mesh...';
progdlg = uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');

T = BoundaryCleanUp(P,T,Csave); % Remove bad boundary elements

[p,t]=fixmesh(P,T);         % Fix mesh

%p = P; t = T;


%--------------------------------------------------------------------------
% store final mesh results
%--------------------------------------------------------------------------
MESH = createMeshStruct(t,p,MESH,PTS,xyzFun); 

end
