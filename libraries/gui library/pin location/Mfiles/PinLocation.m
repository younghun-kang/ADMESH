function PinLocation(varargin)

%------------------------------------------------------------------------------
% Get gui data & axes handle
%------------------------------------------------------------------------------
gui = % guidata(varargin{1}); % Gui data
pH  = gui.ViewAxes;         % Axes handle

%------------------------------------------------------------------------------
% Check if there is anything currently in the plot window
%------------------------------------------------------------------------------
plotItems = get(pH, 'Children'); % Get all graphics handles in plot

% if no graphics handles exist, return
if isempty(plotItems); return; end

% Do we have a mesh in the plot or an edge structure?
meshPatch   = findobj('Tag','Mesh');

if ~isempty(meshPatch) 
    %--------------------------------------------------------------------------
    % Mesh Option
    %--------------------------------------------------------------------------
    
    % Get lat/lon conversion values
    cpplon = gui.MESH.cpplon;
    cpplat = gui.MESH.cpplat;
    
    % Ask user for coordinates
    % gui.sb.setText('Enter Location... ');
    
    msg = 'Check the appropriate coordinate system?';
    choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'Geographic','Cartesian (Meters)'},'DefaultOption',1,'Icon','Warning');
    drawnow; pause(1e-5);

    prompt = {'Enter Latitude/Y:','Enter Longitude/X:'};
    dlg_title = 'ADMESH';
    num_lines = 1;
    def = {'',''};
    input = inputdlg(prompt,dlg_title,num_lines,def);
    
    % User cancels
    if isempty(input); % gui.sb.setText('Ready'); return; end
    
    Y = str2double(input{1});
    X = str2double(input{2});
    
    switch choice,
        case 'Geographic',
            
            if isempty(cpplon)
                
                errordlg('Domain must have been in geographic coordinates.','ADMESH')
                
                drawnow
                return;
                
            end
            
            [X,Y] =  Geo2Meters(X,Y,cpplon,cpplat);
            
    end % switch
    
    % Determine if point is in bounding box
    xLim = get(pH,'xlim'); yLim = get(pH,'ylim');
    
    if X < xLim(1) || X > xLim(2) || Y < yLim(1) || Y > yLim(2)
        
        errordlg('Point out of bounds.','ADMESH')
        
        drawnow
        % gui.sb.setText('Ready');
        return;
        
    end
    
    
    % Plot marker
    h = plot(...
        X,Y,...
        'ko',...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','r',...
        'MarkerSize',10);
    
    hcmenu = uicontextmenu;
    uimenu(hcmenu,'Label','Delete','Callback',{@DeletePinLocation,h});
    set(h,'uicontextmenu',hcmenu)
    uistack(h,'top')
    
    % gui.sb.setText('Ready');
    
else
    %--------------------------------------------------------------------------
    % Edge Structure Option
    %--------------------------------------------------------------------------
    
    % Get lat/lon conversion values
    cpplon = gui.PTS.cpplon;
    cpplat = gui.PTS.cpplat;
    
    axes(pH); hold on % Tell matlab where to plot
    
    % Aske user if they would like to select pin location or enter
    % coordinates
    
    %*****************UNDER CONSTRUCTION********************************
%     choice = questdlg('Select pin location or enter coordinates?', ...
%         'ADMESH', ...
%         'Select Location','Enter Location','Select Location');
%     drawnow; pause(1e-5);
    
    choice = 'Enter Location';

    % Response
    switch choice
        case 'Select Location'
            
            % gui.sb.setText('Select node locations, right-click when complete...');
            
            % Create button down function
            set([gui.ViewAxes ; meshPatch],'Buttondownfcn',{@moveNodes,'select node',[]})
            
            
            
        case 'Enter Location'
            
            % Ask user for coordinates
            % gui.sb.setText('Enter Location...');
            
            msg = 'Check the appropriate coordinate system?';
            choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
                'Options',{'Geographic','Cartesian (Meters)'},'DefaultOption',1,'Icon','Warning');
            drawnow; pause(1e-5);
            
            prompt = {'Enter Latitude/Y:','Enter Longitude/X:'};
            dlg_title = 'ADMESH';
            num_lines = 1;
            def = {'',''};
            input = inputdlg(prompt,dlg_title,num_lines,def);
            
            % User cancels
            if isempty(input); % gui.sb.setText('Ready'); return; end
            
            
            Y = str2double(input{1});
            X = str2double(input{2});

            switch choice,
                case 'Geographic',
                    
                    if isempty(cpplon)
                        
                        errordlg('Domain must have been in geographic coordinates.','ADMESH')
                        
                        drawnow
                        % gui.sb.setText('Ready');
                        return;
                        
                    end
                    
                    [X,Y] =  Geo2Meters(X,Y,cpplon,cpplat);
                    
            end % switch
            
            % Determine if point is in bounding box
            xLim = get(pH,'xlim'); yLim = get(pH,'ylim');
            
            if X < xLim(1) || X > xLim(2) || Y < yLim(1) || Y > yLim(2)
                
                errordlg('Point out of bounds.','ADMESH')
                
                drawnow
                
                % gui.sb.setText('Ready');
                
                return;
                
            end
            
            
            % Plot marker
            h = plot(...
                X,Y,...
                'ko',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',10);
            
            hcmenu = uicontextmenu;
            uimenu(hcmenu,'Label','Delete','Callback',{@DeletePinLocation,h});
            set(h,'uicontextmenu',hcmenu)
            uistack(h,'top')
                        
            % gui.sb.setText('Ready');
            
    end
    
end

%--------------------------------------------------------------------------
% Delete Pin Location subroutine
%--------------------------------------------------------------------------
    function DeletePinLocation(varargin)
        
        delete(varargin{3})
        
    end
end