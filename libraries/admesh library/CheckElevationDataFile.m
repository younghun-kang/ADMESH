function CheckElevationDataFile(app)

filename = app.ElevationDataFilename;
if isempty(filename)
    msg = 'No elevation data is specified.';
    choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');
    % Without this if statement, the subapp is closed immediately
    if strcmp(choice,'OK')
        delete(app);
    end

elseif ~isempty(filename) && exist(filename,'file') == 0
    msg = ['Specified elevation data "',filename, '" cannot be loaded.'];
    choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');
    % Without this if statement, the subapp is closed immediately
    if strcmp(choice,'OK')
        delete(app);
    end

end