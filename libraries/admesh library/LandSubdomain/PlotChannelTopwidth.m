function PlotChannelTopwidth(app,fL,fR,varargin)
hObj = varargin{1};
event = varargin{2};

CP = event.IntersectionPoint(1:2);

H = CP(2);
hp = findobj(hObj,'Tag','ChannelTopWidth');
delete(hp);
hp = plot(hObj,[fL(H),fR(H)],[H H],'b','Tag','ChannelTopWidth');

set(hObj,'UserData',[fL(H),fR(H)]);



end