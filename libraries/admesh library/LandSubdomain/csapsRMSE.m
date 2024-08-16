function varargout = csapsRMSE(x,y,SmoothingRMSE,varargin)
%==========================================================================
% Parser optional arguments
%==========================================================================
%--------------------------------------------------------------------------
% Set default values
%--------------------------------------------------------------------------
xq = [];
id_fixedPoints = [];
%--------------------------------------------------------------------------
% Get optional values
%--------------------------------------------------------------------------
SkipNext = false;
for i = 1 : length(varargin)
    if SkipNext
        SkipNext = false;
        continue;
    end
    ivar = varargin{i};
    if isnumeric(ivar) && isvector(ivar)
        xq = ivar;
    elseif strcmpi(ivar,'fixedpoints')
        id_fixedPoints = find(ismember([x,y],varargin{i+1},'rows'));
        SkipNext = true;
    else
        error('Incorrect use of optional arguments.');
    end
end

%% ========================================================================
w = ones(length(x),1);
w(id_fixedPoints) = 1e8;
    
SP1 = 0;
SP2 = 1;
while 1
    if SP2 - SP1 < 1e-16
        warning(['Iteration cannot reach to the stopping criterion. ',...
            'It is stopped and result with RMSE closest to the given value.']);
        break;
    end
    SmoothingParam = (SP1 + SP2)*0.5;

    Smooth_x = csaps(x,y,SmoothingParam,x,w);
    
    RMSE = sqrt( sum((y(:) - Smooth_x(:)).^2)/length(y));
    if RMSE < SmoothingRMSE
        SP2 = SmoothingParam;
    else
        SP1 = SmoothingParam;
    end
    
    if abs(RMSE - SmoothingRMSE)/SmoothingRMSE < 1e-3
        break;
    end
end

if isempty(xq)
    varargout{1} = csaps(x,y,SmoothingParam);
else
    varargout{1} = csaps(x,y,SmoothingParam,xq);
end
varargout{2} = SmoothingParam;



