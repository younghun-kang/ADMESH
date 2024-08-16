function CreateLandWaterMask(app)

%% ========================================================================
% Check input parameters
%==========================================================================
[status,msg] = CheckParamsExtractChannelsWaterApp(app);
if ~status
    uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');
    return;
end

%% ========================================================================
% Parse input arguments
%==========================================================================
A_min = app.MinAreaEditField.Value;

%% ========================================================================
% Parameter setting
%==========================================================================
%%========================================================================
% Set boundary points from shorelines and set up parameters
%==========================================================================
%--------------------------------------------------------------------------
% Construct polyshape with filtering small islands
%--------------------------------------------------------------------------
progdlg = uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',...
    'Creating Land-Water Mask...','Indeterminate','on');
PTS = app.PTS;
Shoreline1 = struct2table(PTS.Poly,'AsArray',true);
x = Shoreline1.x(:);
y = Shoreline1.y(:);
x = cellfun(@(x) x(~isnan(x)),x,'UniformOutput',0);
y = cellfun(@(x) x(~isnan(x)),y,'UniformOutput',0);
A = cellfun(@(x,y) polyarea(x,y),x,y,'UniformOutput',0);
A = cell2mat(A);
I = A > A_min;
pgon = polyshape(x(I),y(I));
app.pgon_ext = polyshape(x(1),y(1));

app.pgon = subtract(app.pgon_ext,pgon);

PlotLandWaterMask(app);

close(progdlg);

msg = ['Carefully check if blue area correctly indicates water part. ',...
    'If it is reversed, apply "Switch Land-Water Mask".'];
uiconfirm(app.UIFigure,msg,'ADMESH',...
    'Options',{'OK'},'DefaultOption',1,'Icon','Warning');

CheckButtonsExtractChannelsWaterApp(app);

