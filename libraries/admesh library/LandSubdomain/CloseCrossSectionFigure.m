function CloseCrossSectionFigure(h,varargin)
hObj = varargin{1};
event = varargin{2};

delete(h);
delete(hObj);