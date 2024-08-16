function CheckExternalBoundary(app)

PTS = app.PTS;
if isempty(PTS)
    return;
end

X = vertcat(PTS.Poly(2:end).x);
Y = vertcat(PTS.Poly(2:end).y);
X = X(~isnan(X));
Y = Y(~isnan(Y));

IN = insidepoly(X,Y,PTS.Poly(1).x,PTS.Poly(1).y);

flag = all(IN);

if ~flag
    msg = ['Some internal boundaries (green line) are not included in external boundary (black line). ',...
        'Check the input file.'];
    choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');
end
