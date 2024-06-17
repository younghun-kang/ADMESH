function SwitchLandWaterMask(app)

progdlg = uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',...
    'Switching Land-Water Mask...','Indeterminate','on');

app.pgon = subtract(app.pgon_ext,app.pgon);

PlotLandWaterMask(app);
close(progdlg);