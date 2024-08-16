function UpdateSpatialResolutionEditField(app)

switch lower(app.SpatialResolutionDropDown.Value)
    case {'low',1}
        value = 1;
        
    case {'medium',2}
        value = 2;
        
    case {'high',3}
        value = 3;
        
    case {'high (2x)',4}
        value = 4;
        
    case {'high (3x)',5}
        value = 5;
end

app.SpatialResolutionEditField.Value = value;