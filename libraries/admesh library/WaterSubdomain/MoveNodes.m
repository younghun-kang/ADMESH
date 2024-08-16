function MoveNodes(app)

MESH = app.MESH;
t = MESH.ConnectivityList;
p = MESH.Points(:,1:2);

cp = app.UIAxes.CurrentPoint(1,1:2);

id = knnsearch(p,cp);



