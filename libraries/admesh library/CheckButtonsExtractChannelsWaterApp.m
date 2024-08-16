function CheckButtonsExtractChannelsWaterApp(app)

if ~isempty(app.pgon) && ~isempty(app.pgon_ext)
    app.SwitchLandWaterMaskButton.Enable = 'on';
    app.ExtractInternalConstraintsButton.Enable = 'on';
    if ~isempty(app.Constraints)
        app.WriteADMESHinputfileButton.Enable = 'on';
    end
end

