function CheckElevationDataFile(app)

filename = app.ElevationDataFilename;
if isempty(filename) || exist(filename,'file') == 0
    msg = 'No elevation data is loaded.';
    choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');
    % Without this if statement, the subapp is closed immediately
    if strcmp(choice,'OK')
        delete(app);
    end
end