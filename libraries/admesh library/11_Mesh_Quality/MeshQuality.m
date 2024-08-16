function [varargout] = MeshQuality(p,t,Plot,type)
% MeshQuality - 
%
% Triangle: Calculates quality of mesh based off of the ratio
% between the radius of the largest inscribed circle (times two) 
% and the smallest circumscribed
% circle: q = 2(rin/rout).
%
% Quad:
%
% Example: 
% 
% xnodes = [-1 -1  1  1 3  2]';
% ynodes = [-1  1  1 -1 1 -1]';
% t = [1 2 3 4; 4 3 5 6];
% [minq,meanq,q]  = MeshQuality([xnodes ynodes],t,0,'Quad')
% 
%
% Syntax:  [[varargout] = MeshQuality(p,t,Plot)
%
% Inputs:
%    p - nodal points
%    t - connectivity list
%    Plot - 1 to plot, 0 to not plot 
%
% Outputs:
%    min_elem  - min element quality
%    mean_elem - mean element quality
%    q         - element quality
%
% Other m-files required: none
% Subfunctions: meshgrid
% MAT-files required: none

% Author: Colton Conroy, Ethan Kubatko
% The Ohio State University
% email address: conroy.51@osu.edu
% August 2013; Last revision: 08-August-2013

%------------- BEGIN CODE --------------

switch type
    
    case 'Triangle'
        
        % Transpose input vectors if needed
        if size(p,1) == 2
            p = p';
        end
        
        % Read in nodal coordinates
        XNODES = p(:,1);
        YNODES = p(:,2);
        
        % Calculate mesh quality
        
        % Calculate length of sides of element
        a = sqrt((XNODES(t(:,2),1)-XNODES(t(:,1),1)).^2+(YNODES(t(:,2),1)-YNODES(t(:,1),1)).^2);
        b = sqrt((XNODES(t(:,3),1)-XNODES(t(:,2),1)).^2+(YNODES(t(:,3),1)-YNODES(t(:,2),1)).^2);
        c = sqrt((XNODES(t(:,1),1)-XNODES(t(:,3),1)).^2+(YNODES(t(:,1),1)-YNODES(t(:,3),1)).^2);
        
        % mesh quality q for each element
        q = ((b+c-a).*(c+a-b).*(a+b-c))./(a.*b.*c);
        
        % Min and Mean of Element qualities
        min_elem = min(min(q(:,1)));
        mean_elem = mean(q(:,1));
        
        if Plot == 1
            
            % Find elements < .1
            min_loc = find(q(:,1)<.1);
            
            xmin_center =  (XNODES(t(min_loc,1),1)+XNODES(t(min_loc,2),1)+XNODES(t(min_loc,3),1))/3;
            ymin_center =  (YNODES(t(min_loc,1),1)+YNODES(t(min_loc,2),1)+YNODES(t(min_loc,3),1))/3;
            
            figure; hold on
            triplot(t,XNODES,YNODES)
            plot(xmin_center, ymin_center, 'g*','markersize',4)

        end
        
    case 'Quad'
                
        % Transpose input vectors if needed
        if size(p,2) == 2
            p = p';
        end

        xnodes = p(1,:);
        ynodes = p(2,:);
        
        % Place each element in a cell
        ti = num2cell(1:size(t,1))';
        
        % Create face vectors for each cell
        fx = cellfun(@(ti)(xnodes(t(ti,[2 3 4 1])) - xnodes(t(ti,[1 2 3 4]))), ti,'uniformoutput',0);
        fy = cellfun(@(ti)(ynodes(t(ti,[2 3 4 1])) - ynodes(t(ti,[1 2 3 4]))), ti,'uniformoutput',0);
        
        % Compute interior angles for each element
        a = cellfun(@(fx,fy)[fx;fy],fx,fy,'uniformoutput',0); % Faces
        b = cellfun(@(fx,fy)[fx(1,[2 3 4 1]);fy(1,[2 3 4 1])],fx,fy,'uniformoutput',0); % Faces before
        theta = cellfun(@(a,b) acos(dot(a,b)/(norm(a)*norm(b))),a,b,'uniformoutput',0);
        
        % Compute element quality
        q = cellfun(@(theta) prod(1-abs((pi/2 - theta)/(pi/2)),2),theta,'uniformoutput',0);
        q = cell2mat(q);
        
        % Min and Mean of Element qualities
        min_elem = min(min(q(:,1)));
        mean_elem = mean(q(:,1));
        
        if Plot
            
            % Find elements < .1
            min_loc = find(q(:,1)<.1);
            
            xmin_center =  (xnodes(t(min_loc,1),1)+xnodes(t(min_loc,2),1)+xnodes(t(min_loc,3),1)+xnodes(t(min_loc,4),1))/4;
            ymin_center =  (ynodes(t(min_loc,1),1)+ynodes(t(min_loc,2),1)+ynodes(t(min_loc,3),1)+ynodes(t(min_loc,4),1))/4;
            
            figure; hold on
            plot(x(t'),y(t'),'b')
            plot(xmin_center, ymin_center, 'g*','markersize',4)
            
            
        end
        
end

if nargout == 2
    
    varargout{1} = min_elem;

    varargout{2} = mean_elem;
    
elseif nargout == 3
    
    varargout{1} = min_elem;

    varargout{2} = mean_elem;
    
    varargout{3} = q; 
    
end
