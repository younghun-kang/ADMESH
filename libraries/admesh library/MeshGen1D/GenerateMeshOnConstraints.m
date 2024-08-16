function GenerateMeshOnConstraints(X,Y,h,Settings,app)

progdlg = uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',...
    'Generating 1D mesh...');
%--------------------------------------------------------------------------
% Set fixed points
%--------------------------------------------------------------------------
fixedPoints = [];
nfixedP = 0;
PI = app.PI;
for i = 1 : length(PI)
    iPI = PI(i);
    nfixedP = nfixedP + 1;
    fixedPoints(nfixedP,:) = [iPI.x(iPI.p(1)),iPI.y(iPI.p(1))];
    nfixedP = nfixedP + 1;
    fixedPoints(nfixedP,:) = [iPI.x(iPI.p(end)),iPI.y(iPI.p(end))];
end

%----------------------------------------------------------------------
% Generate 1D mesh
%----------------------------------------------------------------------
h2 = @(x,y) interp2(X,Y,h,x,y);
for i = 1 : length(PI)
    Mesh1D(i) = MeshGeneration1D(PI(i),fixedPoints,Settings,h2);
    progdlg.Value = i/length(PI);
end
progdlg = uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',...
    'Applying post processes for 1D mesh...','Indeterminate','on');
[Mesh1D,Constraints] = PostMesh1D(Mesh1D,app.PTS.Constraints,app.MinElementSizeEditField.Value);
app.PTS.Constraints = Constraints;

%----------------------------------------------------------------------
% Update constraints
% *This part needs improvement. Currently assuming that the dummy
% constraints are listed at the last part of PTS.Constraints
%----------------------------------------------------------------------
PTS = app.PTS;
IC = find([PTS.Constraints.num] < 0);
% if length(IC) ~= length(PI)
%     warndlg('Something is wrong.','Error');
%     return;
% end

n = length(IC);
for i = 1 : length(Mesh1D)
    iConstraint = PTS.Constraints(IC(i));
    PTS.Constraints(n+i).num = -iConstraint.num;
    PTS.Constraints(n+i).xy = [Mesh1D(i).X,Mesh1D(i).Y];
    PTS.Constraints(n+i).type = 'line';

    x = iConstraint.xy(:,1);
    y = iConstraint.xy(:,2);
    data = iConstraint.data;
    data_interp = zeros(size(PTS.Constraints(n+i).xy,1),size(data,2));
    for j = 1 : size(iConstraint.data,2)        
        f = scatteredInterpolant(x,y,data(:,j));
        data_interp(:,j) = f(PTS.Constraints(n+i).xy);
    end
    PTS.Constraints(n+i).data = data_interp;
    PTS.Constraints(n+i).Kappa = Mesh1D(i).K;
end
PTS.Constraints(IC) = [];

app.PTS = PTS;

app.DummyConstraint = 0;

PlotEdgeStructure(app,.1);