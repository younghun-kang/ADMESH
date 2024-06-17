function CloseChannelsWaterApp(app)

if ~isempty(app.pgon) && ~app.OutputSaved

    msg = 'Changes made in this subapp has been saved. Do you really want to close this app?';
    choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'Close','Cancel'},'DefaultOption',2,'Icon','warning');

    if strcmp(choice,'Cancel')
        return;
    elseif strcmp(choice,'Close')
        h1 = findobj(app.UIAxes,'tag','Land Water Mask');
        h2 = findobj(app.UIAxes,'tag','Internal Constraints');
        delete(h1);
        delete(h2);

        delete(app);
    end

elseif app.OutputSaved == 1 && ~isempty(app.SavedFile) && exist(app.SavedFile,'file')
    msg = 'A new ADMESH input file is created with extracted constraints. Do you want to open it?';
    choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'Yes','No'},'DefaultOption',1,'Icon','question');

    if strcmp(choice,'Yes')
        mainapp = app.MainApp;
        SavedFile = app.SavedFile;
        delete(app);

        % Read in mat file
        [mainapp,status] = ReadMat(SavedFile,mainapp);
        
        if status == 0
            msg = 'There was an error reading in the file.';
            uiconfirm(mainapp.UIFigure,msg,'ADMESH',...
                'Options',{'OK'},'DefaultOption',1,'Icon','Error');
            return
        end
        
        % Plot Boundary
        PlotEdgeStructure(mainapp,.1);
        
        % Convert to cartesian coordinates if needed.
        [mainapp.PTS,status] = CoordinateConversion(mainapp,mainapp.PTS,'auto');

        if status
            mainapp.xyzFun = CoordinateConversion(mainapp,mainapp.xyzFun,'forward',mainapp.PTS.cpplon,mainapp.PTS.cpplat);
            PlotEdgeStructure(mainapp,.1);
        end

    end
end

delete(app);