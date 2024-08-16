function DisplayMeshInfo(varargin)
% PlotBadElements 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% August 2013; Last revision: 12-August-2013

%--------------------------- BEGIN CODE ---------------------------------------

%------------------------------------------------------------------------------
% Deal inputs
%------------------------------------------------------------------------------
if nargin == 1
    
    % Only input is mesh, redisplay all mesh information
    
    [app] = deal(varargin{:});
    
    p = app.MESH.Points;
    t = app.MESH.ConnectivityList;
    action = 'default';
    
else
    
    % Only update action input
    
    [app,p,t,action] = deal(varargin{:});
    
end

% Get plot window handle
% pH = gui.ViewAxes; axes(pH); 

%------------------------------------------------------------------------------
% Display mesh info
%------------------------------------------------------------------------------
switch action
    
    case 'default'
       
        % Get the number of elements
        nElements = size(t,1);
        
        % Get the number of nodes
        nNodes = size(p,1);
        
        % Get element quality info
        [minEQ, meanEQ, EQ] = MeshQuality(p,t,0,'Triangle');
        
        % Display mesh info
        try
            jf = java.text.DecimalFormat;
        catch
            jf = [];
        end
        if isempty(jf)
            str = {...
                ['Number of Elements: ' num2str(nElements)],...
                '',...
                ['Number of Nodes: ' num2str(nNodes)],...
                '',...
                ['Average Element Quality: ' num2str(meanEQ,'%.2f')],...
                '',...
                ['Minimum Element Quality: ' num2str(minEQ,'%.2f')]};
        else
            str = {...
                ['Number of Elements: ' char(jf.format(nElements))],...
                '',...
                ['Number of Nodes: ' char(jf.format(nNodes))],...
                '',...
                ['Average Element Quality: ' num2str(meanEQ,'%.2f')],...
                '',...
                ['Minimum Element Quality: ' num2str(minEQ,'%.2f')]};
        end
        app.ResultsBox.Value = str;
                   
        % Get current minimum element quality threshold
        MinEQ = app.MinEQ;
        
        % Delete any currently displayed bad elements
        delete(findall(app.UIFigure,'Tag','Bad Element'))
        
        % Find bad elements
        ib = find(EQ < MinEQ);
        
        if ~isempty(ib) % Display bad elements
            
            [x,y] = deal(p(:,1),p(:,2));
            
            patchinfo.xdata          = x(t(ib,:))';
            patchinfo.ydata          = y(t(ib,:))';
            patchinfo.FaceColor      = 'r';
            patchinfo.linewidth      = 1.0;
            h                        = patch(patchinfo,'Parent',app.UIAxes);
            set(h,'Tag','Bad Element')
%             uistack(h,'bottom')
            
        end
        
    case 'update'
        
        % Get the number of elements
        nElements           = size(t,1);
        
        % Get the number of nodes
        nNodes              = size(p,1);
        
        % Get element quality info
        [minEQ, meanEQ, EQ] = MeshQuality(p,t,0,'Triangle');
                
        % Delete any currently displayed bad elements
        delete(findall(app.Window,'Tag','Bad Element'))
        
        % Find bad elements
        ib = find(EQ < app.MinEQ);
        
        if ~isempty(ib) % Display bad elements
            
            [x,y] = deal(p(:,1),p(:,2));
            
            patchinfo.xdata          = x(t(ib,:))';
            patchinfo.ydata          = y(t(ib,:))';
            patchinfo.FaceColor      = 'r';
            patchinfo.linewidth      = 1.0;
            h                        = patch(patchinfo);
            set(h,'Tag','Bad Element')
            uistack(h,'bottom')
            
        end
        
        % Display mesh info
        try
            jf = java.text.DecimalFormat;
        catch
            jf = [];
        end
        if isempty(jf)
            str = {...
                ['Number of Elements: ' num2str(nElements)],...
                '',...
                ['Number of Nodes: ' num2str(nNodes)],...
                '',...
                ['Average Element Quality: ' num2str(meanEQ,'%.2f')],...
                '',...
                ['Minimum Element Quality: ' num2str(minEQ,'%.2f')]};
        else
            str = {...
                ['Number of Elements: ' char(jf.format(nElements))],...
                '',...
                ['Number of Nodes: ' char(jf.format(nNodes))],...
                '',...
                ['Average Element Quality: ' num2str(meanEQ,'%.2f')],...
                '',...
                ['Minimum Element Quality: ' num2str(minEQ,'%.2f')]};
        end
        
        set(app.resultsBox,'string',str)
        
end

























