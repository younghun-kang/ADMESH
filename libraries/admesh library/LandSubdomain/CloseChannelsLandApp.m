function CloseChannelsLandApp(app)

% Update mainapp
app.MainApp.PTS = app.PTS;
app.MainApp.MinDrainageArea = app.MinDrainageAreaEditField.Value;

if app.DEMDropDown.Enable
    msg = 'Which DEM data you want to use?';
    choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'Raw','R2D2','Cancel'},'DefaultOption',3,'Icon','question');
    if strcmp(choice,'Raw')
        app.MainApp.xyzFun = app.xyzFun_old;
    elseif strcmp(choice,'R2D2')
        app.MainApp.xyzFun = app.xyzFun_new;
    else
        return;
    end
    app.DEMDropDown.Value = choice;
    PlotR2D2Result(app);
end

delete(app);