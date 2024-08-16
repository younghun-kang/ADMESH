function PlotInternalConstraints(app)

h1 = findobj(app.UIAxes,'tag','Land Water Mask');
h2 = findobj(app.UIAxes,'tag','Internal Constraints');
delete(h1);
delete(h2);

Constraints = app.Constraints;
pgon = app.pgon;

C3x = [];
C3y = [];
for i = 1 : length(Constraints)
    if abs(Constraints(i).num) == 19
        C3x{end+1} = Constraints(i).xy(:,1);
        C3y{end+1} = Constraints(i).xy(:,2);
    end
end
warnstate = warning;
warning('off');
pgon3 = polyshape(C3x,C3y);
warning(warnstate);

plot(app.UIAxes,pgon,'EdgeColor','k','linewidth',1.5,'FaceColor',[0 0.4470 0.7410],'tag','Internal Constraints');
plot(app.UIAxes,pgon3,'FaceColor',.3*[0 0 1],'linewidth',1.5,'tag','Internal Constraints');
for i = 1 : length(Constraints)
    if abs(Constraints(i).num) == 18
        C = 'b';
    elseif abs(Constraints(i).num) == 17
        C = 'r';
    elseif abs(Constraints(i).num) == 19
        C = 'none';
    end
    plot(app.UIAxes,Constraints(i).xy(:,1),Constraints(i).xy(:,2),'color',C,'tag','Internal Constraints');
end
