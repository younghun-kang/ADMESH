function NewPoints = NaNdlm2struct(XYZ,varargin)

%==========================================================================
ip = inputParser;
%...
v2D    = @(x)validateattributes(x,{'numeric'},...
    {'real','size',[nan 2]});
%...
ip.addParameter('Boundary',[],v2D);
%...
ip.parse(varargin{:});
Boundary      = ip.Results.Boundary;
%==========================================================================

x = XYZ(:,1); y = XYZ(:,2);

if ~isempty(Boundary)
    Bx = Boundary(:,1);
    By = Boundary(:,2);
    
    iNaN = [0; find(isnan(x))];
    
    k = 0;
    for i = 1 : length(iNaN)-1
        id = iNaN(i)+1:iNaN(i+1)-1;
        x1 = x(id);
        y1 = y(id);
        
        [inS,onS] = inpolygon(x1,y1,Bx,By);
        if nnz(inS) > 1
            k = k + 1;
            NewPoints{k} = XYZ(id(inS),:);
        end
    end
    
else
    iNaN = [0; find(isnan(x))];
    
    k = 0;
    for i = 1 : length(iNaN)-1
        id = iNaN(i)+1:iNaN(i+1)-1;
        
        if length(id) > 1
            k = k + 1;
            NewPoints{k} = XYZ(id,:);
        end
    end
    
end


end














