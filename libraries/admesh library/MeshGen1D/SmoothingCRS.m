function SmoothFL = SmoothingCRS(S,BP,k)

%--------------------------------------------------------------------------
% Smooth flowlines with CRS algorithm
%--------------------------------------------------------------------------
ysmooth = smooth(S,S.y,'k',k,'nstribs',true);
xsmooth = smooth(S,S.x,'k',k,'nstribs',true);
csmooth = curvature(S,xsmooth,ysmooth);

%--------------------------------------------------------------------------
% Convert results into vector arrays
%--------------------------------------------------------------------------
[~,~,xsmooth1,ysmooth1] = STREAMobj2XY(S,xsmooth,ysmooth);
id = knnsearch([xsmooth,ysmooth],[xsmooth1,ysmooth1],'k',1);
csmooth1 = csmooth(id);

%--------------------------------------------------------------------------
% Re-arrange the arrays
%--------------------------------------------------------------------------
SmoothFL = NaNdlm2struct([xsmooth1,ysmooth1,csmooth1],'Boundary',BP);