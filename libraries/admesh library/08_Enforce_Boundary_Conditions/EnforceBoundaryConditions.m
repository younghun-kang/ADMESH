function h_ic = EnforceBoundaryConditions(h_ic,X,Y,D,PTS,hmax,hmin,sb)
% EnforceBndCnd - Enforces an element size on h_curve and h_lfs for open
% ocean boundaries if needed.
%
% Syntax:  h_ic = Enforce_Boundary_Conditions(h_ic,X,Y,D,hmax,hmin,delta,guiFig)
%
% Inputs:
%    h_curve - curvature mesh size
%    h_lfs - local feature mesh size
%    X - (nxm) x-coordinates to rectangular grid
%    Y - (n,m) y-coordinates to rectangular grid
%    D - Distance Function
%    IBtype - boundary type, looks for 'Open Ocean'.
%    delta - Grid spacing
%    guiFig - handle that identifies the figure
%
% Outputs:
%    h_curve - adjusted curvature mesh size
%    h_lfs - adjusted local feature mesh size
%
% Other m-files required: none
% Subfunctions: InPolygon,bwdist,Compute_Distance_v3
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% August 2013; Last revision: 08-August-2013

%------------- BEGIN CODE --------------

%------------------------------------------------------------------------------
% Set maximum & minimum element size in mesh
%------------------------------------------------------------------------------
h_ic(h_ic>hmax) = hmax;
h_ic(h_ic<hmin) = hmin;
h_ic(D > hmin)  = hmax;

%--------------------------------------------------------------------------
% Check for constraints
%--------------------------------------------------------------------------
if isempty(PTS.Constraints); return; end

%--------------------------------------------------------------------------
% Find which cells contain open ocean boundaries
%--------------------------------------------------------------------------
% Are there any open ocean boundaries?
ix = find(ismember([PTS.Constraints.num],-1));

if ~isempty(ix);
    
    % Find band hmax for open ocean boundary
    ind = find(abs(D) <= hmax);

    %----------------------------------------------------------------------
    % Enforce Element size to be hmax in these areas
    %----------------------------------------------------------------------
    sb.setText('Enforcing open ocean boundary conditions....')
    
    for k = 1:length(ix)
        
        % Create Polygon Structure
        POLY = create_polygon_structure(PTS.Constraints(ix(k)).xy,hmax);
        
        for i = 1:size(POLY.x,1)
            
            % Find indices that fall within polygon
            IN = PointInPolygon(X(ind),Y(ind),POLY.x(i,:),POLY.y(i,:)); % mex version
            
            %hold on
            %plot(POLY.x(i,:),POLY.y(i,:),'r')
            
            h_ic(ind(IN)) = hmax;
            
            IN = PointInPolygon(X(ind),Y(ind),POLY.Circ.x(i,:),POLY.Circ.y(i,:)); % mex version
            
            h_ic(ind(IN)) = hmax;
            
            %hold on
            %plot(POLY.Circ.x(i,:),POLY.Circ.y(i,:),'r')
            
        end
        
    end
    
end

%--------------------------------------------------------------------------
% Find which cells contain external barrier boundaries
%--------------------------------------------------------------------------
% Are there any external boundaries?
ix = find(ismember([PTS.Constraints.num],[3 13 23]));

if ~isempty(ix);
        
    %----------------------------------------------------------------------
    % Enforce Element size to be hmax in these areas
    %----------------------------------------------------------------------
    sb.setText('Enforcing external barrier conditions....')
    
    for k = 1:length(ix)
        
        % Create Polygon Structure
        POLY = create_polygon_structure(PTS.Constraints(ix(k)).xy);
        
        for i = 1:size(POLY.x,1)
            
            % Find band 
            ind = find(abs(D) <= POLY.L(i));
            
            [IN,~,~] = InPolygon(X(ind),Y(ind),POLY.x(i,:),POLY.y(i,:)); % mex version
            
            h_ic(ind(IN)) = POLY.L(i);
            
            %h_ic(ind(IN)) = hmin;
            
            % Find indices that fall within polygon
            [IN,~,~] = InPolygon(X(ind),Y(ind),POLY.Circ.x(i,:),POLY.Circ.y(i,:)); % mex version
            
            h_ic(ind(IN)) = POLY.L(i);
            
            %h_ic(ind(IN)) = hmin;
            
            %hold on
            %plot(POLY.Circ.x(i,:),POLY.Circ.y(i,:),'r')
            
        end
                
    end
    
end


%------------------------------------------------------------------------------
% Find which cells contain internal barrier boundaries
%------------------------------------------------------------------------------

% Are there any internal boundaries?
ix = find(ismember([PTS.Constraints.num],[4 5 24 25]));

if ~isempty(ix);
    
    %----------------------------------------------------------------------
    % Enforce Element size to be the computed element size in these areas
    %----------------------------------------------------------------------
    sb.setText('Enforcing internal barrier conditions....')
    
    for k = 1:length(ix)
        
        % Create Polygon Structure
        POLY = create_polygon_structure(PTS.Constraints(ix(k)).xy);
        
        for i = 1:size(POLY.x,1)
            
            % Find band
            ind = find(abs(D) <= POLY.L(i));
            
            IN = PointInPolygon(X(ind),Y(ind),POLY.x(i,:),POLY.y(i,:)); % mex version
            
            h_ic(ind(IN)) = POLY.L(i);
            % h_ic(ind(IN)) = hmin;
            
            % Find indices that fall within polygon
            IN = PointInPolygon(X(ind),Y(ind),POLY.Circ.x(i,:),POLY.Circ.y(i,:)); % mex version
            
            h_ic(ind(IN)) = POLY.L(i);
            %h_ic(ind(IN)) = hmin;
            
            %hold on
            %plot(POLY.Circ.x(i,:),POLY.Circ.y(i,:),'r')
            
        end
        
    end
    
end

%--------------------------------------------------------------------------
% Find which cells contain river line boundaries
%--------------------------------------------------------------------------
% Are there any internal boundaries?
ix = find(ismember([PTS.Constraints.num],18));

if ~isempty(ix);
    
    % Find band
    ind = find(abs(D) <= hmin);
    
    %----------------------------------------------------------------------
    % Enforce Element size to be the computed element size in these areas
    %----------------------------------------------------------------------
    sb.setText('Enforcing line boundary conditions...')
    
    for k = 1:length(ix)
        
        % Space points equally based on minimum element size
%         q = SpacePolyEqually(PTS.Constraints(ix(k)).xy,hmin);
        q = PTS.Constraints(ix(k)).xy;  % Younghun
        % Create Polygon Structure
        POLY = create_polygon_structure(q);
        
        for i = 1:size(POLY.x,1)
            ind = find(abs(D) <= POLY.L(i)); % Younghun
            IN = PointInPolygon(X(ind),Y(ind),POLY.x(i,:),POLY.y(i,:)); % mex version
            
            h_ic(ind(IN)) = POLY.L(i);
%             h_ic(ind(IN)) = hmin;
            
            %hold on
            %plot(POLY.x(i,:),POLY.y(i,:),'r')
            
            % Find indices that fall within polygon
            IN = PointInPolygon(X(ind),Y(ind),POLY.Circ.x(i,:),POLY.Circ.y(i,:)); % mex version
            
            h_ic(ind(IN)) = POLY.L(i);
%             h_ic(ind(IN)) = hmin;
            
            %hold on
            %plot(POLY.Circ.x(i,:),POLY.Circ.y(i,:),'r')
            
        end
        
    end
    
end

end