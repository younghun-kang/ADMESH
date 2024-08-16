function SetDummyConstraint(app)

PTS = app.PTS;
if isempty(PTS.Constraints) || ~any([PTS.Constraints.num] < 0)
    return;
end
%----------------------------------------------------------------------
% Setup ADMESH input
%----------------------------------------------------------------------
PTS = app.PTS;
PI = app.PI;

IC = find([PTS.Constraints.num] < 0);
if length(IC) ~= length(PI)
    msg = 'Something is wrong.';
    uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');
    return;
end

n = length(IC);
for i = 1 : length(PI)
    p = PI(i).p;
    x1 = PI(i).x(p);
    y1 = PI(i).y(p);
    k1 = PI(i).k(p);
    
    PTS.Constraints(n+i).num = PTS.Constraints(IC(i)).num;
    PTS.Constraints(n+i).xy = [x1,y1];
    PTS.Constraints(n+i).type = 'line';
    PTS.Constraints(n+i).data = PTS.Constraints(IC(i)).data;
    PTS.Constraints(n+i).Kappa = k1;
end
PTS.Constraints(IC) = [];

app.PTS = PTS;