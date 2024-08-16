function XY = sub2xy(sub,x,y)

sub1 = sub;
if ~iscell(sub1)
    sub1 = {sub1};
end

XY = cellfun(@(i) [x(i(:,2)), y(i(:,1))],sub1,'UniformOutput',0);

if length(XY) == 1 && ~iscell(sub)
    XY = XY{1};
end