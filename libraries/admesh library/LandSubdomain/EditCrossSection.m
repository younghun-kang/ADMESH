function EditCrossSection(app,varargin)

if strcmpi(app.EditCrossSection.State,'off')
    return;
end

hObj = varargin{1};
event = varargin{2};

CP = event.IntersectionPoint(1:2);

h = findobj(app.UIAxes,'tag','Line Constraint');
X = h.XData(:);
Y = h.YData(:);

I = knnsearch([X,Y],CP);
CP = [X(I),Y(I)];

Constraints = app.PTS.Constraints;

for i = 1 : length(Constraints)
    if any(ismember(Constraints(i).xy,CP,'rows'))
        j = find(ismember(Constraints(i).xy,CP,'rows'));
        break;
    end
end

CS = app.CrossSection{i}(j,:);
CS = table2struct(CS);

CoordinateSystem = app.CoordinateSystemDropDown.Value;
MaxChannelWidth = app.MaxChannelWidthEditField.Value;
d = linspace(-MaxChannelWidth/2,MaxChannelWidth/2,1e3);
if strcmpi(CoordinateSystem,'Unprojected (decimal degree)')
    scale2m = 1/km2deg(1e-3);
    d = d/scale2m;
else
    scale2m = 1;
end
d = d(:);

[fL,fR] = InverseCSTopo(CS.fz,MaxChannelWidth,CoordinateSystem);

fig = findobj('Type','Figure','Name','Cross-section');
if isempty(fig)
    fig = figure('Name','Cross-section');
    ax = axes(fig);
else
    ax = findobj(fig.Children,'type','axes');
    if length(ax) > 1
        delete(fig.Children);
        ax = axes(fig);
    end
    cla(ax);
    
    h = findobj(app.UIAxes.Children,'tag','edit-cs');
    delete(h);
end

plot(ax,d*scale2m,CS.fz(d),'k');
hold(ax,'on');
plot(ax,CS.db*scale2m,CS.zb,'-^b','Tag','ChannelTopWidth');
ax.ButtonDownFcn = @(hObj,event)PlotChannelTopwidth(hObj,event);
ax.UserData = [];

h1 = plot(app.UIAxes,CP(1),CP(2),'ok','markersize',10,'tag','edit-cs');
h2 = plot(app.UIAxes,CS.fx(d),CS.fy(d),'k','tag','edit-cs');

fig.CloseRequestFcn = @(hObj,event)CloseCrossSectionFigure([h1,h2],hObj,event);

    function PlotChannelTopwidth(varargin)
        hObj = varargin{1};
        event = varargin{2};
        
        CP = event.IntersectionPoint(1:2);
        
        H = CP(2);
        hp = findobj(hObj,'Tag','ChannelTopWidth');
        delete(hp);
        hp = plot(hObj,[fL(H),fR(H)]*scale2m,[H H],'-^b','Tag','ChannelTopWidth');
        
        set(hObj,'UserData',[fL(H),fR(H)]);
        
        app.CrossSection{i}.db{j} = [fL(H),fR(H)];
        app.CrossSection{i}.zb{j} = [H H];
        
%         app.CrossSection{i}(j,:) = struct2table(CS,'AsArray',1);
         
        PlotBanklines(app);
    end

end

