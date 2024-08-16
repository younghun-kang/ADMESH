function FillMeshSink(app)

MESH = app.MESH;
MESH.Points(:,3) = -MESH.Points(:,3);
MESH = fillSinks(MESH);
MESH.Points(:,3) = -MESH.Points(:,3);
app.MESH = MESH;