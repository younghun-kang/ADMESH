function ImportKMLFile(varargin)
% ImportKML - Reads in and saves a polygon or line constraint
%
% Syntax:  ImportKML(varargin)
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
fig = varargin{1};
gui = % guidata(fig);

% Turn off colormap
SetContourStatus(gui.Window,nan,'off');

%---------------------------------------------------------------------
% Read in KML file
%---------------------------------------------------------------------
% Ask user for file info
[filename,pathname] = uigetfile('*.kml','Select KML file...');

% If user cancels
if filename == 0
    
    % gui.sb.setText('Ready')
    % guidata(fig,gui); % Update GUI data
    
    return
    
end

%------------------------------------------------------------------------------
% Erase any sub meshing stuff
%------------------------------------------------------------------------------
%CancelADMESHSubMesh(fig,nan,'c')

% Open file
fid = fopen([pathname filename],'r');

% gui.sb.setText('Reading in KML file...')

% Read file into cell
C = textscan(fid,'%s');

% Close file id
fclose(fid);

%---------------------------------------------------------------------
% What kind of geometry are we reading in?
%---------------------------------------------------------------------
polygon = ~isempty(find(~cellfun(@isempty,strfind(C{1}, '<Polygon>')), 1));
path    = ~isempty(find(~cellfun(@isempty,strfind(C{1}, '<LineString>')), 1));

% Check if file contains polygons and/or lines
if polygon == 0 && path == 0
    msg = 'To be imported a KML file should contain a polygon or path.';
    uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');
    % gui.sb.setText('Ready')
    return
    
end

% Check if KML file contains both polygon and line
if polygon == 1 && path == 1
    msg = ['To be imported successfully polygons and ' ...
        'paths should be seperate KML files.'];
    uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');
    % gui.sb.setText('Ready')
    return
    
end 

%---------------------------------------------------------------------
% If a mesh exists in the program, ask the user if they would like to clear
% it
%---------------------------------------------------------------------
if ~isempty(gui.MESH)
    
    msg = ['To import a KML file ADMESH+ needs to clear the current mesh.',...
        ' What would you like to do?'];
    choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'Save mesh and clear','Clear mesh','Cancel'},'DefaultOption',3,'Icon','Warning');
    switch choice
        
        case 'Save mesh and clear'
            
            % Call save mesh routine
            try
                SaveMesh14File(gui.Window,gui);
            catch err
            
                disp(err) 
                % gui.sb.setText('Ready')
                return
                
            end
            
            % Get plot handles
            plotItems = get(gui.ViewAxes, 'Children');
            
            % Delete from plot
            if ~isempty(plotItems); delete(plotItems); end
            
            % delete color bar if it exists
            colorbar('delete') 
            
            % Remove current legend if it exist
            delete(findall(fig,'Tag','Mesh Info'))
            
            gui.MESH = []; % Mesh data
            
            % guidata(fig,gui)
            
            if ~isempty(gui.PTS)
                PlotEdgeStructure(gui,.1);
            end
            
        case 'Clear mesh'
            
            % Get plot handles
            plotItems = get(gui.ViewAxes, 'Children');
            
            % Delete from plot
            if ~isempty(plotItems); delete(plotItems); end
            
            % delete color bar if it exists
            colorbar('delete')
            
            % Remove current legend if it exist
            delete(findall(fig,'Tag','Mesh Info'))
            
            gui.MESH = []; % Mesh data
            
            % guidata(fig,gui)
            
            if ~isempty(gui.PTS)
                PlotEdgeStructure(gui,.1);
            end

        case 'Cancel'
            
            % gui.sb.setText('Ready')
            return
            
    end
    
end
    
%---------------------------------------------------------------------
% Polygon Case
%---------------------------------------------------------------------
if polygon
    
    % Check for existing PTS data structure
    if ~isempty(gui.PTS)
        
        % Ask user how we're going to import the data
        msg = ['Are you creating a new edge structure or '...
            'appending polygon(s) to the current edge structure file?'];
        choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
            'Options',{'New File','Appending','Cancel'},'DefaultOption',2,'Icon','Warning');
    else
        
        choice = 'New File';
        
    end
    
    % Case
    switch choice,
        
        case 'New File'
            
            gui.FilePath   = [];
            gui.PTS        = [];
            gui.xyzFun     = [];
            gui.MESH       = [];
                        
            % Search for beginning coordinate flags
            idx = find(~cellfun(@isempty,strfind(C{1}, '<coordinates>')))+1;
            
            % Search for ending coordinate flags
            idy = find(~cellfun(@isempty,strfind(C{1}, '</coordinates>')))-1;
            
            % Initialize cell for holding segments
            S = cell(numel(idx),1);
            
            % Initialize area vector
            A = zeros(numel(S),1);
            
            % Store path segments in cell
            gui.sb.Progressbar.setVisible(true)
            set(gui.sb.Progressbar,'Minimum',1,'Maximum',numel(idx),'Value',1)
            % gui.sb.setText('Compiling polygon(s)...')
            
            for k = 1:numel(idx)
                
                % Extract segment
                cellSeg = cellfun(@str2num,C{1}(idx(k):idy(k)),'uniformoutput',0);
                xyz = [cellSeg{:}];
                
                % Convert to cartesian coordinates
                if k == 1
                    [x,y,gui.PTS.cpplon,gui.PTS.cpplat] = Geo2Meters(xyz(1:3:end)',xyz(2:3:end)');
                else
                    [x,y] = Geo2Meters(xyz(1:3:end)',xyz(2:3:end)',gui.PTS.cpplon,gui.PTS.cpplat);
                end
                
                % Make sure point set is unique and closed
                S{k} = [unique([x y],'rows','stable'); [x(1) y(1)]];
                
                % Compute polygon area
                A(k) = polyarea(S{k}(:,1),S{k}(:,2));
                
                set(gui.sb.Progressbar,'Value',k)
                
            end
            
            % Sorting polygons based on area
            [~,i] = sort(A,'descend');
            
            % Create edge structure
            set(gui.sb.Progressbar,'Minimum',1,'Maximum',numel(S),'Value',1)
            for k = 1:numel(S)
                
                gui.PTS.Poly(k).x = S{i(k)}(:,1);
                gui.PTS.Poly(k).y = S{i(k)}(:,2);
                
                set(gui.sb.Progressbar,'Value',k)
                
            end
            
            % guidata(fig,gui)
            
            % Plot Edge Structure
            PlotEdgeStructure(gui,.05);
            
            gui.sb.Progressbar.setVisible(false)
            % gui.sb.setText('Ready')

        case 'Appending'
            
            
            % Search for beginning coordinate flags
            idx = find(~cellfun(@isempty,strfind(C{1}, '<coordinates>')))+1;
            
            % Search for ending coordinate flags
            idy = find(~cellfun(@isempty,strfind(C{1}, '</coordinates>')))-1;
            
            % Initialize cell for holding segments
            S = cell(numel(idx),1);
            
            % initialize warning flaf
            wrnmsg = 0;
            
            % Store path segments in cell
            gui.sb.Progressbar.setVisible(true)
            set(gui.sb.Progressbar,'Minimum',1,'Maximum',numel(idx),'Value',1)
            % gui.sb.setText('Checking path compatability...')
            
            for k = 1:numel(idx)
                
                % Extract section
                xyz =  cell2mat(cellfun(@str2num,C{1}(idx(k):idy(k)),'uniformoutput',0));
                
                % Convert to cartesian coordinates
                [x,y] = Geo2Meters(xyz(:,1),xyz(:,2),gui.PTS.cpplon,gui.PTS.cpplat);
                                
                % Check if segment is in domain
                IN = PointsInDomain2(x,y,gui.PTS);
                
                if all(IN)
                    
                    % Store in cell if all segments are in domain
                    S{k} =  [x y];
                    
                else
                    
                    % Display warning message
                    wrnmsg = 1;
                    
                end
                
                set(gui.sb.Progressbar,'Value',k)
                
            end
            
            gui.sb.Progressbar.setVisible(false)
            
            % Remove empty cells in S
            S(cellfun(@isempty,S)) = [];
            
            % Check if any segments were stored
            if isempty(S)
                
                errordlg(['None of the path segments imported fell '...
                    'with in the edge structure currently loaded.'],'Error');
                % gui.sb.setText('Ready')
                return
                
            end
            
            % Warning message
            if wrnmsg
                
                msg = ['Some paths were not stored because not '...
                    'all points were with in the domain.'];
                uiconfirm(app.UIFigure,msg,'ADMESH',...
                    'Options',{'OK'},'DefaultOption',1,'Icon','Error');
                
            end
            
            msg = 'What type of polygon is this?';
            choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
                'Options',{'Internal Barrier','Internal Boundary','Internal Ridge Line'},'DefaultOption',3,'Icon','Warning');
            % Append to PTS data structure
            gui.sb.Progressbar.setVisible(true)
            set(gui.sb.Progressbar,'Minimum',1,'Maximum',numel(S),'Value',1)
            % gui.sb.setText('Appending to edge structure...')
            
            switch choice
                
                case 'Internal Boundary'
                    
                    % Define constext menu
                    hcmenu = uicontextmenu;
                    uimenu(hcmenu,'Label','Delete Polygon' ,'Callback',@DeletePolygon);
                    
                    for k = 1:numel(S)
                        
                        % Make sure polygon point list is unique & closed
                        xy = [unique(S{k},'rows','stable') ; S{k}(1,:)];
                        
                        gui.PTS.Poly(end+1).x = xy(:,1);
                        
                        gui.PTS.Poly(end).y = xy(:,2);
                        
                        h = plot(gui.PTS.Poly(end).x , gui.PTS.Poly(end).y,'g','LineWidth',LineWidth);
                        set(h,'tag','Edge Structure','UserData',{'Internal Boundary',length(gui.PTS.Poly)})
                        set(h,'uicontextmenu',hcmenu)
                        
                        set(gui.sb.Progressbar,'Value',k)
                        
                    end
                    
                    gui.sb.Progressbar.setVisible(false)
                    % guidata(fig,gui)
                    
                    % gui.sb.setText('Ready')
                    
                case 'Internal Ridge Line'

                    hcmenu = uicontextmenu;
                    uimenu(hcmenu,'Label','Remove BC','Callback',@RemoveBC);
                    
                    
                    for k = 1:numel(S)
                        
                        % Make sure polygon point list is unique & closed
                        xy = [unique(S{k},'rows','stable') ; S{k}(1,:)];

                        gui.PTS.Constraints(end+1).xy = [xy(:,1),xy(:,2)];
                        gui.PTS.Constraints(end).num = 19;
                        gui.PTS.Constraints(end).type = 'Ridge';
                        gui.PTS.Constraints(end).data = [];
                                                
                        h = plot(xy(:,1),xy(:,2),'r');
                        set(h,'tag','Ridge Constraint','UserData',length(gui.PTS.Constraints),'uicontextmenu',hcmenu)
                        uistack(h, 'top')
                        
                        set(gui.sb.Progressbar,'Value',k)
                        
                    end
                    
                    gui.sb.Progressbar.setVisible(false)
                    % guidata(fig,gui)
                    
                    % gui.sb.setText('Ready')
                    
                    
                    
                    
            end
            
        case 'Cancel'
            
            % gui.sb.setText('Ready')
            return
            
    end % switch
    
    % In-case no constraints were added
    if ~isfield(gui.PTS,'Constraints')
        gui.PTS.Constraints = [];
    end
    
    % guidata(fig,gui)
    
    % gui.sb.setText('Ready')
end

%---------------------------------------------------------------------
% Line case
%---------------------------------------------------------------------
if path
    
    % Check for edge structure
    if isempty(gui.PTS)
        errordlg(['To import a path you must have an edge structure '...
            'loaded in ADMESH to append the path to.'],'Error');
        % gui.sb.setText('Ready')
        return
        
    end
    
    % Check for cpplon & cpplat
    if ~isfield(gui.PTS,'cpplon')
        errordlg(['To import a KML file the coordinates of the edge '...
            'structure should have been imported into ADMESH in ' ...
            'geographic format.'],'Error');
        % gui.sb.setText('Ready')
        return
    end
    
    % Check for constraints field
    if ~isfield(gui.PTS,'Constraints')
        gui.PTS.Constraints = []; % Add field if missing
    end
    
    % Search for beginning coordinate flags
    idx = find(~cellfun(@isempty,strfind(C{1}, '<coordinates>')))+1;
    
    % Search for ending coordinate flags
    idy = find(~cellfun(@isempty,strfind(C{1}, '</coordinates>')))-1;
    
    % Initialize cell for holding segments
    S = cell(numel(idx),1);
    
    % Get unit conversion constants
    cpplon = gui.PTS.cpplon;
    cpplat = gui.PTS.cpplat;
    
    % Update status bar
    gui.sb.Progressbar.setVisible(true)
    set(gui.sb.Progressbar,'Minimum',1,'Maximum',numel(idx),'Value',1)
    % gui.sb.setText('Checking path compatability...')
    
    % Store path segments in cell
    for k = 1:numel(idx)
        
        % Extract segment
        cellSeg = cellfun(@str2num,C{1}(idx(k):idy(k)),'uniformoutput',0);
        xyz = [cellSeg{:}];
        
        % Convert to cartesian coordinates
        [x,y] = Geo2Meters(xyz(1:3:end)',xyz(2:3:end)',cpplon,cpplat);
        
        % Check if segment is in domain
        IN = PointsInDomain2(x,y,gui.PTS);
        
        if all(IN)
            
            % Store in cell if all segments are in domain
            S{k} =  [x y];
            
        elseif any(IN)
            
            % Remove coordinates not in domain, replace with NaN
            x(~IN) = nan; y(~IN) = nan;
            
            % Assign to matrix, A
            A = [x,y];
            
            % Split into cells
            pos = [true, isnan(A(:, 1)).', true];
            ini = strfind(pos, [true, false]);
            fin = strfind(pos, [false, true]) - 1;
            C   = cell(length(ini),1);
            for iC = 1:length(ini)
                C{iC} = A(ini(iC):fin(iC), :);
            end
            
            % Store segments that fall in domain, appen to S
            S = [S; C]; %#ok<AGROW>
        end
        
        set(gui.sb.Progressbar,'Value',k)
        
    end
    
    gui.sb.Progressbar.setVisible(false)
    
    % Remove empty cells in S
    S(cellfun(@isempty,S)) = [];
    
    % Check if any segments were stored
    if isempty(S)
        errordlg(['None of the path segments imported fell '...
            'with in the edge structure currently loaded.'],'Error');
        % gui.sb.setText('Ready')
        return
    end
    
    
    % Distance tolerance
    tol = .0001;
    
    % Segment channel lines
    gui.sb.Progressbar.setVisible(true)
    set(gui.sb.Progressbar,'Minimum',1,'Maximum',numel(S),'Value',1)
    % gui.sb.setText('Processing path segments...')
    
    for k = 1:numel(S)
        
        % Store end points
        EP1 = [S{k}(1,1)  , S{k}(1,2)];
        EP2 = [S{k}(end,1), S{k}(end,2)];
        
        % Indices to test
        j = 1:numel(S); j(j == k) = [];
        
        for i = 1:numel(j)
            
            % Test first end point.
            ix = find( ...
                sqrt((S{j(i)}(:,1) - EP1(1)).^2 + (S{j(i)}(:,2) - EP1(2)).^2)...
                <= tol ,1,'first' );
            
            % If we don't find anything move on
            if ~isempty(ix)
                
                % Make matching points equal
                S{j(i)}(ix(1),:) = EP1;
                
                % Break up line segment
                S{j(i)} = [S{j(i)}(1:ix(1),:); nan nan; S{j(i)}(ix(1):end,:)];
                
            end
            
            % Test second end point.
            ix = find( ...
                sqrt((S{j(i)}(:,1) - EP2(1)).^2 + (S{j(i)}(:,2) - EP2(2)).^2)...
                <= tol ,1,'first' );
            
            % If we don't find anything move on
            if ~isempty(ix)
                
                % Make matching points equal
                S{j(i)}(ix(1),:) = EP2;
                
                % Break up line segment
                S{j(i)} = [S{j(i)}(1:ix(1),:); nan nan; S{j(i)}(ix(1):end,:)];
                
            end
            
        end
        
        set(gui.sb.Progressbar,'Value',k)
        
    end
    
    set(gui.sb.Progressbar,'Minimum',1,'Maximum',numel(S),'Value',1)
    % gui.sb.setText('Appending to edge structure...')
    
    % Define the context menu items and install their callbacks for constraints
    hcmenu = uicontextmenu;
    uimenu(hcmenu,'Label','Remove BC','Callback',@RemoveBC);
    
    % Append to PTS data structure
    for k = 1:numel(S)
        
        xy = unique(S{k},'rows','stable');
        
        if any(isnan(xy(:,1)))
            
            % Break up array into segments from NaN seperation
            i1 = all(~isnan(xy),2); i2 = i1(:)';
            idx = [strfind([~i2(1),i2],[0 1]); strfind([i2, ~i2(end)],[1 0])];
            seg = mat2cell(xy(i1,:),diff(idx)+1,size(xy,2));
            
            % Store segments
            for j = 1:length(seg)
                
                gui.PTS.Constraints(end+1).xy = seg{j};
                gui.PTS.Constraints(end).num = 18;
                gui.PTS.Constraints(end).type = 'Line';
                gui.PTS.Constraints(end).data = [];
                
                h = plot(seg{j}(:,1),seg{j}(:,2),'r');
                set(h,'tag','Line Constraint','UserData',length(gui.PTS.Constraints),'uicontextmenu',hcmenu)
                uistack(h, 'top')
                
                % Add to sub mesh structure
                if isfield(gui,'subPTS')
                    
                    gui.subPTS.Constraints(end+1).xy = seg{j};
                    gui.subPTS.Constraints(end).num = 18;
                    gui.subPTS.Constraints(end).type = 'Line';
                    gui.subPTS.Constraints(end).data = [];

                end
                
            end
            
        else
            
            gui.PTS.Constraints(end+1).xy = [xy(:,1),xy(:,2)];
            gui.PTS.Constraints(end).num = 18;
            gui.PTS.Constraints(end).type = 'Line';
            gui.PTS.Constraints(end).data = [];
            
            % Add to sub mesh structure
            if isfield(gui,'subPTS')
                
                gui.subPTS.Constraints(end+1).xy = [xy(:,1),xy(:,2)];
                gui.subPTS.Constraints(end).num = 18;
                gui.subPTS.Constraints(end).type = 'Line';
                gui.subPTS.Constraints(end).data = [];
                
            end
            
            h = plot(xy(:,1),xy(:,2),'r');
            set(h,'tag','Line Constraint','UserData',length(gui.PTS.Constraints),'uicontextmenu',hcmenu)
            uistack(h, 'top')
            
        end
        
        set(gui.sb.Progressbar,'Value',k)
        
    end
    
    gui.sb.Progressbar.setVisible(false)
    
    % guidata(fig,gui)
    
    % Plot Edge Structure
    %PlotEdgeStructure(gui,.1);
    
    % gui.sb.setText('Ready')

end

end