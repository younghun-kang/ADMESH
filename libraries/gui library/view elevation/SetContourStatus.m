function SetContourStatus(app,varargin)
% ADMESh gui function for viewing elevation data
%
% writting by Dustin West (west.425@osu.edu)

%--------------------------------------------------------------------------
% Get GUI data
%--------------------------------------------------------------------------
if length(varargin) == 1
    status = varargin{1};
    app.ContoursDropDown.Value = status;
else
    status = app.ContoursDropDown.Value;
end

%--------------------------------------------------------------------------
% Check for data
%--------------------------------------------------------------------------
if ~strcmpi(status,'off')
    if strcmpi(status,'bathy/topo')

        h1 = findobj(app.UIAxes,'tag','Edge Structure');

        if ~isempty(h1) && isempty(app.xyzFun)
            msg = 'No vertical elevation data are identified.';
            uiconfirm(app.UIFigure,msg,'ADMESH',...
                'Options',{'OK'},'DefaultOption',1,'Icon','Error');
            SetContourStatus(app,'Off');
            return;
        end
        
    end

    if strcmpi(status,'Manning''s n') && ~isfield(app.MESH,'Attributes')
        msg = 'No Manning''n n value data are identified.';
        uiconfirm(app.UIFigure,msg,'ADMESH',...
            'Options',{'OK'},'DefaultOption',1,'Icon','Error');
        SetContourStatus(app,'Off');
        return;
    end
end

%--------------------------------------------------------------------------
% Are we turning on or off?
%--------------------------------------------------------------------------
switch lower(status)
    
    case {'bathy/topo','manning''s n'}
        
        % Enable dropdown menus
        app.ColormapDropDown.Enable = 'On';

        % Get current colormap
        cmap_name = app.ColormapDropDown.Value;
                        
        % Are we displaying elevation nodes or triangulated elevation
        h1 = findobj(app.UIAxes,'tag','Edge Structure');
        h2 = findobj(app.UIAxes,'tag','Mesh');
        
        msg = 'Setting colormap...';
        progdlg = uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');
        
        if ishandle(h1) == 1 % Plot elevation nodes
            
            Field = -app.xyzFun.Values;

            % Get colormap
            if strcmpi(cmap_name,'land & sea')

                cmap = demcmap(Field,512,flipud(seacolor(256)),landcolor(256));
                
            elseif strcmpi(cmap_name,'jet')

                cmap = flipud(jet(256));
                
            elseif strcmpi(cmap_name,'parula')

                cmap = flipud(parulacolor(256));
                
            end
            
            msg = 'Preparing to display...';
            progdlg = uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');
            
            % Plot
            if isa(app.xyzFun,'scatteredInterpolant')
                
                
                delete(findobj(app.UIAxes,'tag','elevation plot'))
                
                % Get x & y limits
                xmin = min(vertcat(app.PTS.Poly(:).x));
                xmax = max(vertcat(app.PTS.Poly(:).x));
                ymin = min(vertcat(app.PTS.Poly(:).y));
                ymax = max(vertcat(app.PTS.Poly(:).y));
                
                % Compute image offset
                offset = max((xmax-xmin),(ymax-ymin));
                
                % Re-compute limits with offset
                xmin = xmin - offset*.05;
                xmax = xmax + offset*.05;
                ymin = ymin - offset*.05;
                ymax = ymax + offset*.05;
                
                % Find indices that fall within limits
                ixy = ...
                    app.xyzFun.Points(:,1) >= xmin & ...
                    app.xyzFun.Points(:,1) <= xmax & ...
                    app.xyzFun.Points(:,2) >= ymin & ...
                    app.xyzFun.Points(:,2) <= ymax;
                
                % Assign vectors
                x = app.xyzFun.Points(ixy,1)';
                y = app.xyzFun.Points(ixy,2)';
                z = Field(ixy)';
                
                % Designate colormap
                colormap(app.UIAxes,cmap)
                
                % Plot points
                eH = patch(...
                    'parent',app.UIAxes,...
                    'xdata',x,...
                    'ydata',y,...
                    'cdata',z,...
                    'facecolor','none',...
                    'edgecolor','none',...
                    'marker','sq',...
                    'markerfacecolor','flat',...
                    'tag','elevation plot');
                
                uistack(eH,'bottom')

                drawnow;
                
            else
                
                delete(findobj(app.UIAxes,'tag','elevation plot'))
                delete(findobj(app.UIAxes,'tag','colorbar'))
                
               % Get x & y limits
                xmin = min(vertcat(app.PTS.Poly(:).x));
                xmax = max(vertcat(app.PTS.Poly(:).x));
                ymin = min(vertcat(app.PTS.Poly(:).y));
                ymax = max(vertcat(app.PTS.Poly(:).y));
                
                % Compute image offset
                offset = max((xmax-xmin),(ymax-ymin));
                
                % Re-compute limits with offset
                xmin = xmin - offset*.05;
                xmax = xmax + offset*.05;
                ymin = ymin - offset*.05;
                ymax = ymax + offset*.05;
                
                % Assign vectors
                x = app.xyzFun.GridVectors{1};
                y = app.xyzFun.GridVectors{2};
                
                % Find indices that do not fall within limits
                ix = ~(x >= xmin & x <= xmax);
                iy = ~(y >= ymin & y <= ymax);
                
                % Remove points that fall outside range
                Field(ix,:) = [];
                Field(:,iy) = [];
                
                % Designate colormap
                
                
                eH = imagesc(app.UIAxes,x(~ix),y(~iy),Field','tag','elevation plot');
                colormap(app.UIAxes,cmap);
                axis(app.UIAxes,'xy');
                
                % Put elevation plot at the bottom
                for i = 1 : length(app.UIAxes.Children)
                    if strcmpi(app.UIAxes.Children(i).Tag,'elevation plot')
                        app.UIAxes.Children = app.UIAxes.Children([1:i-1,i+1:end,i]);
                        break;
                    end
                end
                
            end
            
            drawnow;
            close(progdlg);
            
        elseif ishandle(h2) == 1
            
            switch lower(status)
                case 'bathy/topo'
                    Field = -app.MESH.Points(:,3);
                case 'manning''s n'
                    Field = app.MESH.Attributes.ManningN;
            end
            

            % Get colormap
            if strcmpi(cmap_name,'land & sea')

                cmap = demcmap(Field,512,seacolor(256),landcolor(256));
                
            elseif strcmpi(cmap_name,'jet')

                cmap = flipud(jet(256));
                
            elseif strcmpi(cmap_name,'parula')

                cmap = flipud(parulacolor(256));
                
            end
            
            msg = 'Preparing to display...';
            progdlg = uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');
            
            delete(h2);

            verts = app.MESH.Points;                        % vertices
            faces = app.MESH.ConnectivityList;              % connectivity list
            
            trisurf(faces,verts(:,1),verts(:,2),Field,...
                'Parent',app.UIAxes,...
                'Tag','Mesh',...
                'EdgeColor','none',...
                'FaceColor','interp');

            colormap(app.UIAxes,cmap);

        end

        colorbar(app.UIAxes,'Tag','Mesh');
        drawnow;

    case 'off'
        
        % Check box on
        app.ColormapDropDown.Enable = 'off';
        
        % Are we displaying elevation nodes or triangulated elevation
        h1 = findobj(app.UIAxes,'tag','Edge Structure');
        h2 = findobj(app.UIAxes,'tag','Mesh');
        
        if ishandle(h1) == 1
            
            delete(findobj(app.UIAxes,'tag','elevation plot'))
            delete(findobj(app.UIAxes,'tag','colorbar'))
            
        elseif ishandle(h2) == 1
            
            delete(h2);
            if ~isempty(app.MESH)
                PlotMesh(app,.1);
            end
        end

        colorbar(app.UIAxes,'off');
        drawnow;

end

    function y = landcolor(n)
        %LANDCOLOR Land colormap
        %
        %	Author: Francois Beauducel <beauducel@ipgp.fr>
        %	$Revision: 1.0.0 $   $Date: 2012/05/17 11:22:44 $
        
        J = [ ...
            0.095678 0.53427 0.21682
            0.15785 0.5979 0.23274
            0.21286 0.64673 0.2514
            0.26411 0.68789 0.27268
            0.32959 0.72416 0.31308
            0.39794 0.75695 0.36038
            0.46153 0.7871 0.40624
            0.52108 0.81516 0.45135
            0.57702 0.84152 0.49547
            0.62973 0.86645 0.53891
            0.67946 0.89016 0.58187
            0.72647 0.91282 0.62427
            0.77095 0.93455 0.66619
            0.81306 0.95546 0.70772
            0.85292 0.97563 0.7489
            0.89066 0.99514 0.78976
            0.88379 0.98595 0.77038
            0.86389 0.96758 0.73236
            0.84615 0.94972 0.69623
            0.8303 0.93233 0.66186
            0.81612 0.91536 0.6291
            0.80341 0.8988 0.59784
            0.79201 0.8826 0.56795
            0.78191 0.86676 0.53946
            0.7729 0.85123 0.51224
            0.76479 0.83602 0.48615
            0.75747 0.8211 0.46111
            0.75084 0.80645 0.43704
            0.74506 0.79206 0.41414
            0.73981 0.77792 0.39211
            0.73501 0.76401 0.37089
            0.73068 0.75033 0.35052
            0.72683 0.73685 0.33106
            0.72042 0.72074 0.31228
            0.71032 0.70085 0.29417
            0.69761 0.67821 0.27694
            0.68489 0.65558 0.26026
            0.67235 0.63313 0.24418
            0.65997 0.61082 0.22889
            0.64775 0.58874 0.21406
            0.63568 0.56689 0.19983
            0.62376 0.54527 0.18622
            0.61197 0.52391 0.17299
            0.60033 0.50283 0.16046
            0.58881 0.48203 0.14832
            0.57742 0.46151 0.13667
            0.56616 0.44133 0.12555
            0.55502 0.4214 0.11472
            0.54398 0.4019 0.10456
            0.53306 0.38266 0.094633
            0.52226 0.36382 0.085242
            0.51155 0.3453 0.076179
            0.50095 0.32714 0.067515
            0.49045 0.30938 0.059259
            0.48005 0.29193 0.051294
            0.46973 0.27495 0.043796
            0.45951 0.25823 0.0365
            0.44938 0.24206 0.029715
            0.43934 0.22609 0.023063
            0.42938 0.21074 0.016949
            0.41951 0.19556 0.010917
            0.40971 0.18105 0.0054326
            0.4 0.16667 0
            ];
                
        l = length(J);
        if nargin < 1
            n = 256;
        end
        y = interp1(1:l,J,linspace(1,l,n),'*linear');
        
    end

    function y = parulacolor(n)
        % parulaCOLOR parula colormap
        
        J = [ ...
            0.2081    0.1663    0.5292
            0.2116    0.1898    0.5777
            0.2123    0.2138    0.6270
            0.2081    0.2386    0.6771
            0.1959    0.2645    0.7279
            0.1707    0.2919    0.7792
            0.1253    0.3242    0.8303
            0.0591    0.3598    0.8683
            0.0117    0.3875    0.8820
            0.0060    0.4086    0.8828
            0.0165    0.4266    0.8786
            0.0329    0.4430    0.8720
            0.0498    0.4586    0.8641
            0.0629    0.4737    0.8554
            0.0723    0.4887    0.8467
            0.0779    0.5040    0.8384
            0.0793    0.5200    0.8312
            0.0749    0.5375    0.8263
            0.0641    0.5570    0.8240
            0.0488    0.5772    0.8228
            0.0343    0.5966    0.8199
            0.0265    0.6137    0.8135
            0.0239    0.6287    0.8038
            0.0231    0.6418    0.7913
            0.0228    0.6535    0.7768
            0.0267    0.6642    0.7607
            0.0384    0.6743    0.7436
            0.0590    0.6838    0.7254
            0.0843    0.6928    0.7062
            0.1133    0.7015    0.6859
            0.1453    0.7098    0.6646
            0.1801    0.7177    0.6424
            0.2178    0.7250    0.6193
            0.2586    0.7317    0.5954
            0.3022    0.7376    0.5712
            0.3482    0.7424    0.5473
            0.3953    0.7459    0.5244
            0.4420    0.7481    0.5033
            0.4871    0.7491    0.4840
            0.5300    0.7491    0.4661
            0.5709    0.7485    0.4494
            0.6099    0.7473    0.4337
            0.6473    0.7456    0.4188
            0.6834    0.7435    0.4044
            0.7184    0.7411    0.3905
            0.7525    0.7384    0.3768
            0.7858    0.7356    0.3633
            0.8185    0.7327    0.3498
            0.8507    0.7299    0.3360
            0.8824    0.7274    0.3217
            0.9139    0.7258    0.3063
            0.9450    0.7261    0.2886
            0.9739    0.7314    0.2666
            0.9938    0.7455    0.2403
            0.9990    0.7653    0.2164
            0.9955    0.7861    0.1967
            0.9880    0.8066    0.1794
            0.9789    0.8271    0.1633
            0.9697    0.8481    0.1475
            0.9626    0.8705    0.1309
            0.9589    0.8949    0.1132
            0.9598    0.9218    0.0948
            0.9661    0.9514    0.0755
            0.9763    0.9831    0.0538
            ];
        
        l = length(J);
        if nargin < 1
            n = 256;
        end
        y = interp1(1:l,J,linspace(1,l,n),'*linear');
        
        
    end

    function y = seacolor(n)
        %SEACOLOR Sea colormap adapted from NGDC ETOPO1
        %
        %	Author: Francois Beauducel <beauducel@ipgp.fr>
        
        J = [ ...
            0.0392         0    0.4745
            0.1020         0    0.5373
            0.1020         0    0.5373
            0.1490         0    0.5961
            0.1490         0    0.5961
            0.1059    0.0118    0.6510
            0.1059    0.0118    0.6510
            0.0627    0.0235    0.7059
            0.0627    0.0235    0.7059
            0.0196    0.0353    0.7569
            0.0196    0.0353    0.7569
            0    0.0549    0.7961
            0    0.0549    0.7961
            0    0.0863    0.8235
            0    0.0863    0.8235
            0    0.1176    0.8471
            0    0.1176    0.8471
            0    0.1529    0.8745
            0    0.1529    0.8745
            0.0471    0.2667    0.9059
            0.0471    0.2667    0.9059
            0.1020    0.4000    0.9412
            0.1020    0.4000    0.9412
            0.0745    0.4588    0.9569
            0.0745    0.4588    0.9569
            0.0549    0.5216    0.9765
            0.0549    0.5216    0.9765
            0.0824    0.6196    0.9882
            0.0824    0.6196    0.9882
            0.1176    0.6980    1.0000
            0.1176    0.6980    1.0000
            0.1686    0.7294    1.0000
            0.1686    0.7294    1.0000
            0.2157    0.7569    1.0000
            0.2157    0.7569    1.0000
            0.2549    0.7843    1.0000
            0.2549    0.7843    1.0000
            0.3098    0.8235    1.0000
            0.3098    0.8235    1.0000
            0.3686    0.8745    1.0000
            0.3686    0.8745    1.0000
            0.5412    0.8902    1.0000
            0.5412    0.8902    1.0000
            0.7373    0.9020    1.0000
            ];
        
        l = length(J);
        if nargin < 1
            n = 256;
        end
        y = interp1(1:l,J,linspace(1,l,n),'*linear');
        
        
    end

end
