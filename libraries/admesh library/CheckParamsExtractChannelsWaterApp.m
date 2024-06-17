function [status,msg] = CheckParamsExtractChannelsWaterApp(app)

%% ========================================================================
% Parse input arguments
%==========================================================================
A_min = app.MinAreaEditField.Value;
dx = app.BackgroundGridResolmEditField.Value;
delta_cw = app.MaxChannelWidthmEditField.Value;

%% ========================================================================
% Check input parameters
%==========================================================================
status = 1;
msg = {'Issues with the following input parameters:',''};
if A_min < 0
    msg{end+1} = '- "Min Area" must be a non-negative value.';
end
if dx <= 0 || delta_cw < 0
    msg{end+1} = '- "Background Grid Resolution" must be a positive value.';
end
if delta_cw < 0
    msg{end+1} = '- "Max Channel Width" must be a non-negative value.';
end
if length(msg) > 2
    status = 0;
end

