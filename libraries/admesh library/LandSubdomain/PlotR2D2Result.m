function PlotR2D2Result(app)

switch app.DEMDropDown.Value
    case 'R2D2'
        app.MainApp.xyzFun = app.xyzFun_new;
    case 'Raw'
        app.MainApp.xyzFun = app.xyzFun_old;
end

SetContourStatus(app.MainApp,'Bathy/Topo');

