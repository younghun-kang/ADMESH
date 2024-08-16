function PlotBanklines(app)

CS = app.CrossSection;

x = [];
y = [];
for i = 1 : length(CS)
    D = vertcat(CS{i}.db);
    I = cellfun(@(x) length(x) < 2,D);
    D(I) = {[0 0]};
    
    iCS_fx = CS{i}.fx;
    iCS_fy = CS{i}.fy;
    
    x = [x; cellfun(@(f,x) f(x(1)),iCS_fx(:),D); nan];
    y = [y; cellfun(@(f,x) f(x(1)),iCS_fy(:),D); nan];
        
    x = [x; cellfun(@(f,x) f(x(2)),iCS_fx(:),D); nan];
    y = [y; cellfun(@(f,x) f(x(2)),iCS_fy(:),D); nan];
    
end

h = findobj(app.UIAxes,'tag','banklines');
delete(h);

c = [150, 75, 0]/255;
h = plot(app.UIAxes,x,y,'Color',c);
set(h,'tag','banklines');

