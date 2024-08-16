function Compute1DProjection(app)

PTS = app.PTS;
if isempty(PTS.Constraints) || ~any([PTS.Constraints.num] < 0)
    app.PI = [];
    app.DummyConstraint = 0;
    return;
end

IC = find([PTS.Constraints.num] < 0);
Points = {PTS.Constraints(IC).xy};
SmoothingRMSE = app.SmoothingRMSEEditField.Value;

%----------------------------------------------------------------------
% Add fixed points for junctions
%----------------------------------------------------------------------
fixedPoints = [];
nfixedP = 0;
for j = 1 : length(Points)
    nfixedP = nfixedP + 1;
    fixedPoints(nfixedP,:) = [Points{j}(1,1),Points{j}(1,2)];
    nfixedP = nfixedP + 1;
    fixedPoints(nfixedP,:) = [Points{j}(end,1),Points{j}(end,2)];
end

clear PI;
k = 0;
for i = 1 : length(Points)
    x = Points{i}(:,1);
    y = Points{i}(:,2);

    iPI = ComputePathCurvature([x,y],fixedPoints,SmoothingRMSE);
    if isempty(iPI)
        PI(i).p = 0;
    else
        PI(i) = iPI;
    end
end

hmin = app.MinElementSizeEditField.Value;
p = {PI.p};
I = cellfun(@(x) x(end) == 0,p);
PI(I) = [];
PTS.Constraints(IC(I)) = [];

app.PTS = PTS;
app.PI = PI;

if app.DummyConstraint == 0
    app.DummyConstraint = 1;
end

end