function Z = CreateElevationGrid(X,Y,xyzFun,Settings,UIFigure)


%------------------------------------------------------------------------------
% Check for On/Off & Available data
%------------------------------------------------------------------------------
OnOff = strcmp(Settings.B.Status,'On') || strcmp(Settings.T.Status,'On');

if OnOff && ~isempty(xyzFun) % If true, interpolate bathymetry to grid
    
    msg = 'Interpolating bathymetry onto background mesh...';
    uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');

    %---------------------------------------------------------------------------
    % Interpolate bathymetry (xyzFun) onto background mesh (X,Y)
    %---------------------------------------------------------------------------   
    Z = xyzFun(X',Y')';
        
    % Interpolate NaN's if they exist
    if any(isnan(Z(:)))
        Z = inpaint_nans(Z);
    end

else
    
    Z = [];

end


end



