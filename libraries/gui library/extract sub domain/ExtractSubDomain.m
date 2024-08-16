function ExtractSubDomain(varargin)
% ExtractSubDomain - Extracts subdomain from PTS data structure
%
% Syntax:  ExtractSubDomain(varargin)
%
% Inputs:
%    None
%
% Outputs:
%    None
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% Last revision: 29-April-2014

%---------------------------------------------------------------------
% Begin Code
%---------------------------------------------------------------------

%---------------------------------------------------------------------
% Get GUI Data
%---------------------------------------------------------------------
fig = varargin{1}; gui = % guidata(fig);

%---------------------------------------------------------------------
% Check if any data isloaded
%---------------------------------------------------------------------
if isempty(gui.PTS) && isempty(gui.MESH)
    errordlg('No file is currently loaded.','Error');
    return
end

%---------------------------------------------------------------------
% Check for mesh and no edge structure
%---------------------------------------------------------------------
if isempty(gui.PTS) && ~isempty(gui.MESH)
    
    % Ask user what to do
    msg = [...
        'ADMESH needs to import your current mesh into an edge structure '...
        'file. Select Continue to proceed. Select '...
        'Cancel to quit.'];
    choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'Continue','Cancel'},'DefaultOption',1,'Icon','Warning');
    
    drawnow; pause(.005);
    
    % To continue or not to continue.
    switch choice
        
        case 'Continue' % Create edge structure file
            
            ImportMesh(gui);
        
            gui = % guidata(gui.Window);
            
        case 'Cancel'
            return
    end
    
end

% Assign variables
PTS     = gui.PTS;
xyzFun  = gui.xyzFun;

% Turn off colormap
SetContourStatus(gui.Window,nan,'off');

%---------------------------------------------------------------------
% Check for mesh in plot window
%---------------------------------------------------------------------
pH = gui.ViewAxes;           % Plot Window handle
h = findobj(pH,'Tag','Mesh');   % Mesh plot tag

% Check for mesh in plot window
if ~isempty(h)

    % Ask user if they are sure they want to delete the mesh
    msg = ['ADMESH needs to clear this mesh. '...
        'Are you sure you want to clear this mesh?'];
    choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'Clear','Cancel'},'DefaultOption',2,'Icon','Warning');
    drawnow; pause(0.05);  % this innocent line prevents the Matlab hang
    
    if strcmp(choice,'Cancel'),% gui.sb.setText('Ready');return,end

    plotItems = get(pH, 'Children');
    
    if ~isempty(plotItems); delete(plotItems); end
    
    colorbar('delete') % Delete current colorbar
    
    % Remove current legend if it exist
    delete(findall(fig,'Tag','Mesh Info'))
    
    gui.MESH = []; % Mesh data
    gui.hmin = []; % Minimum element size
    
    % guidata(fig,gui)
    
    % Clear mesh from screen & Plot Boundary
    PlotEdgeStructure(PTS,.1);
    
end

%---------------------------------------------------------------------
% Ask user to draw polygon
%---------------------------------------------------------------------
bpoly = DrawBoundingPolygon(gui.Window,pH,'b');

if isempty(bpoly)
    % gui.sb.setText('Ready')
    return
end

%---------------------------------------------------------------------
% Extract sub domain
%---------------------------------------------------------------------
% gui.sb.setText('Extracting sub domain')

nSeg = numel(PTS.Poly);

subD = cell(nSeg,1); % Cell for holding sub domain

ntc  = false(nSeg,1); % Flag indicating whether polygon needs to be completed

gui.sb.ProgressBar.setVisible(true)
set(gui.sb.ProgressBar, 'Minimum',1, 'Maximum',nSeg, 'Value',1)

for k = 1:nSeg
    
    % Find points in bpoly
    IN = PointInPolygon(PTS.Poly(k).x(1:end-1),PTS.Poly(k).y(1:end-1),bpoly(:,1),bpoly(:,2));
    
    % If any points are in
    if any(IN)
        
        % Assign polygon points
        x = PTS.Poly(k).x(1:end-1);
        y = PTS.Poly(k).y(1:end-1);
        
        % Is segment complete?
        if any(~IN) % No
            
            % Create a column vector representing indices in x & y
            ix = (1:length(x))';
            
            % Set indices that do not below to zero
            ix(IN == 0) = 0;
            
            % Shift indices so we can extract one column vector of points
            shiftN = abs(find(ix == 0,1,'last') - length(ix));
            
            % Shift
            ix = circshift(ix,shiftN);
            
            % Remove zeros
            ix(ix == 0) = [];
            
            % Save to subD
            subD{k} = [x(ix),y(ix)];
            
            % Set flag
            ntc(k) = true;
            
        else % Yes
            
            % Save to subD
            subD{k} = [x,y];
            
        end

    end
    
    set(gui.sb.ProgressBar,'Value',k)

end

gui.sb.ProgressBar.setVisible(false)

% Remove empty cells
emptyCells = cellfun(@isempty,subD);

subD(emptyCells) = [];
ntc(emptyCells)  = [];

%---------------------------------------------------------------------
% Create tPTS data structure
%---------------------------------------------------------------------
tPTS.Poly = repmat(struct('x',1,'y',1), length(ntc)-1, 1 );

for i = 1:length(ntc)
    
    xy = unique(subD{i},'rows','stable');
    
    if ntc(i) % If we still need to complete
        tPTS.Poly(i).x = xy(:,1);
        tPTS.Poly(i).y = xy(:,2);
    else
        tPTS.Poly(i).x = xy(:,1);
        tPTS.Poly(i).y = xy(:,2);
        
        % Close polygon
        tPTS.Poly(i).x(end+1) = tPTS.Poly(i).x(1);
        tPTS.Poly(i).y(end+1) = tPTS.Poly(i).y(1);
    end
    
end

% Find polygons that need completed
ntc = find(ntc);

%--------------------------------------------------------------------------
% Complete polygons
%--------------------------------------------------------------------------
for i = 1:length(ntc)
    
    if ntc(i) == 1 % Add open ocean boundary
                
        while 1
                        
            %----------------------------------------------------------------------
            % Set axis prpoerties
            %----------------------------------------------------------------------
            PlotEdgeStructure(tPTS,gui,gui.sb,pH,.1)
            
            zoomH = findobj('Tag','Zoom Tag');
            panH  = findobj('Tag','Pan Tag');
            
            % Set all other toolbar tools off
            set([zoomH panH],'State','off')
            
            % Have user select open ocean boundary
            ttPTS = DrawOpenOceanBoundary(tPTS,gui,pH,gui.sb);
            
            % Plot new PTS
            PlotEdgeStructure(ttPTS,gui,gui.sb,pH,.1)
            
            % gui.sb.setText('Continue?')
            
            % Ask the user if they like their selection
            msg = [...
                'Select Continue to keep the coordinates of the '...
                'boundary you entered. Select Re-Try to enter different '...
                'coordinates.'];
            choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
                'Options',{'Continue','Re-Try','Quit'},'DefaultOption',1,'Icon','Warning');
            
            drawnow; pause(0.05);  % this innocent line prevents the Matlab hang
            
            switch choice
                
                case 'Continue'
                    
                    tPTS = ttPTS; clear ttPTS 
                    break
                    
                case 'Re-Try'
                    
                    clear ttPTS
                    
                    % Clear window
                    plotItems = get(pH, 'Children');
                    
                    if ~isempty(plotItems); delete(plotItems); end
                    
                    PlotEdgeStructure(tPTS,gui,gui.sb,pH,.1)
                    
                case 'Quit'
                    
                    clear ttPTS tPTS
                    
                    % Clear window
                    plotItems = get(pH, 'Children');
                    
                    if ~isempty(plotItems); delete(plotItems); end
                    
                    PlotEdgeStructure(gui,.1);
                    
                    % gui.sb.setText('Ready')
                    
                    return
                    
            end
            
            
        end
        
    else
        
        while 1
            
            %----------------------------------------------------------------------
            % Set axis prpoerties
            %----------------------------------------------------------------------
            PlotEdgeStructure(tPTS,.1)
            
            zoomH = findobj('Tag','Zoom Tag');
            panH  = findobj('Tag','Pan Tag');
            
            % Set all other toolbar tools off
            set([zoomH panH],'State','off')
            
            % Have user select open ocean boundary
            ttPTS = CloseInternalBoundary(tPTS,fig,pH,ntc(i));
            
            % Plot new PTS
            PlotEdgeStructure(ttPTS,.1)
            
            % gui.sb.setText('Continue?')
            
            % Ask the user if they like their selection
            msg = [...
                'Select Continue to keep the coordinates of the '...
                'boundary you entered. Select Re-Try to enter different '...
                'coordinates.'];
            choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
                'Options',{'Continue','Re-Try','Quit'},'DefaultOption',1,'Icon','Warning');
            drawnow; pause(0.05);  % this innocent line prevents the Matlab hang
            
            switch choice
                
                case 'Continue'
                    
                    tPTS = ttPTS; clear ttPTS 
                    break
                    
                case 'Re-Try'
                    
                    clear ttPTS
                    
                    % Clear window
                    plotItems = get(pH, 'Children');
                    
                    if ~isempty(plotItems); delete(plotItems); end
                    
                    PlotEdgeStructure(tPTS,.1)
                    
                case 'Quit'
                    
                    clear ttPTS tPTS
                    
                    % Clear window
                    plotItems = get(pH, 'Children');
                    
                    if ~isempty(plotItems); delete(plotItems); end
                    
                    PlotEdgeStructure(PTS,.1);
                    
                    % gui.sb.setText('Ready')
                    
                    return
                    
            end
            
            
        end
        
    end
    
end

% Save cpplon & cpplat if they exist
if isfield(PTS,'cpplon')
    
   tPTS.cpplon = PTS.cpplon;
   tPTS.cpplat = PTS.cpplat;
   
else
    
   tPTS.cpplon = [];
   tPTS.cpplat = [];
   
end

% Assign PTS, clear temporary data structure
PTS = tPTS; clear tPTS

%--------------------------------------------------------------------------
% Check for islands outside main domain
%--------------------------------------------------------------------------
% gui.sb.setText('Checking domain for exterior islands...')
gui.sb.ProgressBar.setVisible(true)
nSegments = numel(PTS.Poly);
set(gui.sb.ProgressBar, 'Minimum',1, 'Maximum',nSegments, 'Value',1)

ir = false(1,nSegments); % indices to remove

for i = 2:nSegments
    
    set(gui.sb.ProgressBar,'Value',i)
        
    IN = PointInPolygon(PTS.Poly(i).x, PTS.Poly(i).y,PTS.Poly(1).x, PTS.Poly(1).y);
    
    if any(~IN); ir(i) = true; end 
    
end

gui.sb.ProgressBar.setVisible(false)

PTS.Poly(ir) = []; % Remove exterior islands

% gui.sb.setText( 'Save File...')

%--------------------------------------------------------------------------
% Save as a .mat file
%--------------------------------------------------------------------------
[filename,pathname] = uiputfile('*.mat','Save File As');

if filename ~=0

    gui.FilePath = [pathname filename];% Save filename/pathname
    
    % gui.sb.setText( 'Saving File...')
    
    save([pathname, filename], 'PTS', 'xyzFun','-v7.3');
    
    gui.FilePath = [pathname filename];
    
    % Update GUI workspace
    gui.PTS     = PTS;
    gui.xyzFun  = xyzFun;
    
    % gui.sb.setText( 'Plotting in cartesian coordinates')
    
    % guidata(gui.Window,gui); 
    
    % Plot Boundary
    PlotEdgeStructure(gui,.1);
    
    % gui.sb.setText( 'Ready')
    
else
    
    % Clear window
    plotItems = get(pH, 'Children');
    
    if ~isempty(plotItems); delete(plotItems); end
    
    set(pH ,...
        'XLim',[0 1],...
        'YLim',[0 1])
    drawnow
    
    % Update GUI workspace
    gui.PTS     = [];
    gui.xyzFun  = [];
    
    % gui.sb.setText( 'Ready')
    
    % guidata(gui.Window,gui);
    
end

end