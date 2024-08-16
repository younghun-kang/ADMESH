function D2MA = DistanceToMA(MA,MaskCell,dd_ID,dx,UIFigure)

msg = 'Computing distance to Medial Axis...';
progdlg = uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');

%% ========================================================================
% Separate 1D & 2D area based on local feature size
%==========================================================================

if ~isempty(MA.BranchNodes)
id_MA = vertcat(MA.BranchNodes{:});
if islogical(id_MA)
    id_MA = find(id_MA);
end
else
    id_MA = [];
end

if ~iscell(MaskCell)
    MaskCell = {MaskCell};
end
%--------------------------------------------------------------------------
% More accurate but too slow..
%--------------------------------------------------------------------------
% d2MA = bwdistgeodesic(MASK,id_MA);
% d2MA(isinf(d2MA)) = 0;
% d2MA = d2MA*dx;

%--------------------------------------------------------------------------
% Use domain decomposition to speed up
%--------------------------------------------------------------------------
id_MA = vertcat(MA.BranchNodes{:});
[id_MA1,id_MA2] = ind2sub(MA.Size,id_MA);

D2MA = sparse(MA.Size(1),MA.Size(2));
for iDD = 1 : size(dd_ID,1)
    I = dd_ID{iDD,1};
    J = dd_ID{iDD,2};
    
    width = max(J) - min(J) + 1;
    height = max(I) - min(I) + 1;
    width = ceil(width/2);
    height = ceil(height/2);
    
    
    I = min(I) - height : max(I) + height;
    J = min(J) - width : max(J) + width;
    
    I(I < 1 | I > MA.Size(1)) = [];
    J(J < 1 | J > MA.Size(2)) = [];
    
    K = id_MA1 >= min(I) & id_MA1 <= max(I) & id_MA2 >= min(J) & id_MA2 <= max(J);
    
    I2 = id_MA1(K);
    J2 = id_MA2(K);
    
    I2 = I2 - min(I) + 1;
    J2 = J2 - min(J) + 1;
    
    MaskMA1 = false(length(I),length(J));
    K2 = sub2ind(size(MaskMA1),I2,J2);
    MaskMA1(K2) = 1;
    
    D2MA1 = full(D2MA(I,J));
    for iMask = 1 : length(MaskCell)
        MASK = full(MaskCell{iMask}(I,J));
        d2MA = bwdistgeodesic(MASK,MaskMA1);
        d2MA(isnan(d2MA)) = 0;
        %     d2MA(isinf(d2MA)) = 0;
        d2MA = d2MA*dx;
        
        K = abs(D2MA1) == 0 | (abs(d2MA) > 0 & abs(D2MA1) > abs(d2MA));
        
        D2MA1(K) = d2MA(K);
    end
    
    D2MA(I,J) = D2MA1;
    
    progdlg.Indeterminate = 'off';
    progdlg.Value = iDD/size(dd_ID,1);
end

