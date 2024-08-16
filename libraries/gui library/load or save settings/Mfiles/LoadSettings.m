function LoadSettings(Settings,app)

if isempty(Settings); return; end

%------------------------------------------------------------------------------
% Set Resolution Setting
%------------------------------------------------------------------------------
switch lower(Settings.Res)
    case {'low',1}
        Res = 'Low';
        
    case {'medium',2}
        Res = 'Medium';
        
    case {'high',3}
        Res = 'High';
        
    case {'high 2x','high (2x)',4}
        Res = 'High (2x)';
        
    case {'high 3x','high (3x)',5}
        Res = 'High (3x)';
        
end
app.SpatialResolutionDropDown.Value = Res;
UpdateSpatialResolutionEditField(app);

%------------------------------------------------------------------------------
% Load Element Sizes
%------------------------------------------------------------------------------
app.MaxElementSizeEditField.Value = Settings.hmax;
app.MinElementSizeEditField.Value = Settings.hmin;

%------------------------------------------------------------------------------
% Load Curvature Settings
%------------------------------------------------------------------------------
switch lower(Settings.K.Status)
    case 'on'
        app.BoundaryCurvatureDropDown.Value = 'On';
        app.BoundaryCurvatureEditField.Editable = 'On';
        app.BoundaryCurvatureEditField.Value = Settings.K.Value;
        
    case 'off'
        app.BoundaryCurvatureDropDown.Value = 'Off';
        app.BoundaryCurvatureEditField.Editable = 'Off';
end

%------------------------------------------------------------------------------
% Load LFS Settings
%------------------------------------------------------------------------------
switch lower(Settings.R.Status)
    case 'on'
        app.LocalFeatureSizeDropDown.Value = 'On';
        app.LocalFeatureSizeEditField.Value = Settings.R.Value;

    case 'off'
        app.LocalFeatureSizeDropDown.Value = 'Off';
        app.LocalFeatureSizeEditField.Editable = 'Off';
end

%------------------------------------------------------------------------------
% Load Bathymetry Settings
%------------------------------------------------------------------------------
switch lower(Settings.B.Status)
    case 'on'
        app.ElevationGradientsDropDown.Value = 'On';
        app.ElevationGradientsEditField.Value = Settings.B.Value;
       
    case 'off'
        app.ElevationGradientsDropDown.Value = 'Off';
        app.ElevationGradientsEditField.Editable = 'Off';
end

%------------------------------------------------------------------------------
% Load Tidal Settings
%------------------------------------------------------------------------------
switch lower(Settings.T.Status)
    case 'on'
        app.DominateTideDropDown.Value = 'On';
        app.DominateTideEditField.Value = Settings.T.Value;
       
    case 'off'
       app.DominateTideDropDown.Value = 'Off';
       app.DominateTideEditField.Editable = 'Off';
        
end

%------------------------------------------------------------------------------
% Load Grading Settings
%------------------------------------------------------------------------------
switch lower(Settings.G.Status)
    case 'on'
        app.MeshGradingDropDown.Value = 'On';
        app.MeshGradingEditField.Editable = 'On';
        app.MeshGradingEditField.Value = Settings.G.Value;
                
    case 'off'
        app.MeshGradingDropDown.Value = 'Off';
        app.MeshGradingEditField.Editable = 'Off';

end

%------------------------------------------------------------------------------
% Load View Settings
%------------------------------------------------------------------------------
switch lower(Settings.View.Status)
    case 'on'
        app.ViewMeshGenerationDropDown.Value = 'On';

    case 'off'
        app.ViewMeshGenerationDropDown.Value = 'Off';

end

%--------------
if isfield(Settings,'DummyConstraint')
    app.DummyConstraint = Settings.DummyConstraint;
end

if isfield(Settings,'ElevationDataFilename')
    app.ElevationDataFilename = Settings.ElevationDataFilename;
end

if isfield(Settings,'SmoothingRMSE')
    app.SmoothingRMSEEditField.Value = Settings.SmoothingRMSE;
end

if isfield(Settings,'MinDrainageArea')
    app.MinDrainageArea = Settings.MinDrainageArea;
end

end