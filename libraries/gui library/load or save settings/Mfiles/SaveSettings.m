function [Settings,AdvSettings] = SaveSettings(app)

%------------------------------------------------------------------------------
% Save Resolution Setting
%------------------------------------------------------------------------------
Settings.Res = app.SpatialResolutionDropDown.Value;

%------------------------------------------------------------------------------
% Save Element Sizes
%------------------------------------------------------------------------------
Settings.hmax = app.MaxElementSizeEditField.Value;
Settings.hmin = app.MinElementSizeEditField.Value;
%------------------------------------------------------------------------------
% Save Curvature Settings
%------------------------------------------------------------------------------
switch lower(app.BoundaryCurvatureDropDown.Value)
    case 'on'
        Settings.K.Status = 'On';
        Settings.K.Value  = app.BoundaryCurvatureEditField.Value;
    case 'off'
        Settings.K.Status = 'Off';
        Settings.K.Value  = nan;
end

%------------------------------------------------------------------------------
% Save LFS Settings
%------------------------------------------------------------------------------
switch lower(app.LocalFeatureSizeDropDown.Value)
    case 'on'
        Settings.R.Status = 'On';
        Settings.R.Value  = app.LocalFeatureSizeEditField.Value;
    case 'off'
        Settings.R.Status = 'Off';
        Settings.R.Value  = nan;
end

%------------------------------------------------------------------------------
% Save Bathymetry Settings
%------------------------------------------------------------------------------
switch lower(app.ElevationGradientsDropDown.Value)
    case 'on'
        Settings.B.Status = 'On';
        Settings.B.Value  = app.ElevationGradientsEditField.Value;
    case 'off'
        Settings.B.Status = 'Off';
        Settings.B.Value  = nan;
end

%------------------------------------------------------------------------------
% Save Tidal Settings
%------------------------------------------------------------------------------
switch lower(app.DominateTideDropDown.Value)
    case 'off'
        Settings.T.Status   = 'Off';
        Settings.T.Value    = nan;
    case {'m2','s2','n2','k2'}
        Settings.T.Status   = 'On';
        Settings.T.Value    = app.DominateTideEditField.Value;
        Settings.T.type     = app.DominateTideDropDown.Value;
end

%------------------------------------------------------------------------------
% Save Grading Settings
%------------------------------------------------------------------------------
switch lower(app.MeshGradingDropDown.Value)
    case 'on'
        Settings.G.Status = 'On';
        Settings.G.Value  = app.MeshGradingEditField.Value;
    case 'off'
        Settings.G.Status = 'Off';
        Settings.G.Value  = nan;
end

%------------------------------------------------------------------------------
% Save Viewing Settings
%------------------------------------------------------------------------------
switch lower(app.ViewMeshGenerationDropDown.Value)
    case 'on'
        Settings.View.Status = 'On';
    case 'off'
        Settings.View.Status = 'Off';
end

Settings.DummyConstraint = app.DummyConstraint;
Settings.SmoothingRMSE = app.SmoothingRMSEEditField.Value;
Settings.MinDrainageArea = app.MinDrainageArea;

%------------------------------------------------------------------------------
% Save Advanced Settings
%------------------------------------------------------------------------------
AdvSettings.ttol       = app.ttolEditField.Value;
AdvSettings.Fscale     = app.FscalEditField.Value;
AdvSettings.deltat     = app.deltatEditField.Value;
AdvSettings.geps_scale = app.geps_scaleEditField.Value;
AdvSettings.niter      = app.niterEditField.Value;
AdvSettings.qold       = app.qoldEditField.Value;



end