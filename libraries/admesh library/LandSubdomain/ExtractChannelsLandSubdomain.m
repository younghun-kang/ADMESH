function ExtractChannelsLandSubdomain(app,event)

% First call
if strcmpi(event,'firstcall')
    filename = app.ElevationDataFilename;
    if isempty(filename) || isempty(app.xyzFun)
        msg = 'No elevation data is read.';
        uiconfirm(app.UIFigure,msg,'ADMESH',...
            'Options',{'OK'},'DefaultOption',1,'Icon','Error');
        return;
    end
    [~,~,ext] = fileparts(filename);
    
    if ~any(strcmpi(ext,{'.tif','.tiff'}))
        msg = 'The elevation data should be .tif or .tiff format.';
        uiconfirm(app.UIFigure,msg,'ADMESH',...
            'Options',{'OK'},'DefaultOption',1,'Icon','Error');
        return;
    end
    
    PTS = app.PTS;
    if ~isempty(PTS.Constraints) && any([PTS.Constraints(:).num] == 18)
        msg = 'Open channel constraints are found. Do you want to replace them?';
        choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
            'Options',{'Yes','No'},'DefaultOption',2,'Icon','Warning');
        
        if strcmpi(choice,'no')
            return;
        end
        
        I = [PTS.Constraints(:).num] == 18;
        PTS.Constraints(I) = [];
        app.PTS = PTS;
        
        PlotEdgeStructure(app,.1);
    end
    
    progdlg = uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',...
        'Extracting open channels from DEM...','Indeterminate','on');
    [FD,A] = ExtractOpenChannels(filename);
    app.TTC.FD = FD;
    app.TTC.A = A;
end

if isempty(app.TTC)
    % Do nothing and return if not called with button first
    return;
end
    
FD = app.TTC.FD;
A = app.TTC.A;

progdlg = uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',...
    'Extracting open channels from DEM...','Indeterminate','on');
warnStruct = warning;
warning('off');
FL1 = klargestconncomps(STREAMobj(FD,A > app.MinDrainageAreaEditField.Value));
warning(warnStruct);

if isempty(FL1)
    msg = ['There are only 0 connected components in the stream network. ',...
        'Use smaller threshold to get stream network.'];
    uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');
    
    app.TTC.FL = {};
    
    h = findobj(app.UIAxes,'tag','landchannel');
    delete(h);
    return;
end

BP = [app.PTS.Poly.x(:), app.PTS.Poly.y(:)];
[sx,sy] = STREAMobj2XY(FL1);
if ~isempty(app.PTS.cpplon) && ~isempty(app.PTS.cpplat)
    % Convert to XY (meters)
    [sx,sy] = Geo2Cart(sx,sy,app.PTS.cpplon,app.PTS.cpplat);
end
FL = NaNdlm2struct([sx,sy],'Boundary',BP);

% Remove existing constraints extracted by TopoToolbox
PTS = app.PTS;
if ~isempty(PTS.Constraints) && isfield(PTS.Constraints(1),'type')
    type = {PTS.Constraints.type};
    I = cellfun(@(x) strcmp(x,'TopoToolbox Channel'),type);
    PTS.Constraints(I) = [];
end
% Add new channels to constraints
n = length(PTS.Constraints);
for i = 1 : length(FL)
    PTS.Constraints(n+i).num = -18;
    PTS.Constraints(n+i).xy = FL{i};
    PTS.Constraints(n+i).type = 'TopoToolbox Channel';
    PTS.Constraints(n+i).data = [];
    PTS.Constraints(n+i).Kappa = [];
end
app.PTS = PTS;

% Visualize
progdlg = uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',...
    'Draw extracted open channels...','Indeterminate','on');
FL2 = cellfun(@(x) [x; nan(1,2)],FL,'UniformOutput',0);
FL2 = vertcat(FL2{:});

h = findobj(app.UIAxes,'tag','landchannel');
delete(h);
h = plot(app.UIAxes,FL2(:,1),FL2(:,2),'b');
drawnow;

set(h,'tag','landchannel');

% Update UI
app.NumR2D2LinesEditField.Enable = 1;
app.RunR2D2SmoothingButton.Enable = 1;

close(progdlg);
