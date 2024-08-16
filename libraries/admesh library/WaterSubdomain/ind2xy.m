function XY = ind2xy(Ind,Size,x,y)

Ind1 = Ind;
if ~iscell(Ind1)
    Ind1 = {Ind1};
end
[I,J] = cellfun(@(x) ind2sub(Size,x),Ind1,'UniformOutput',0);
XY = cellfun(@(i,j) [x(j(:)), y(i(:))],I,J,'UniformOutput',0);

if length(XY) == 1 && ~iscell(Ind)
    XY = XY{1};
end