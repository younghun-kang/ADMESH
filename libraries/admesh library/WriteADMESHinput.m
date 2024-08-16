function WriteADMESHinput(app)

Constraints = app.Constraints;

progdlg = uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',...
    'Select a file...','Indeterminate','on');
% Ask the user to select a file
% f_dummy = warndlg('Select a file...','ADMESH'); %create a dummy figure so that uigetfile doesn't minimize our GUI
f_dummy = figure('Position',[-100 -100 0 0]);
[filename, pathname] = uiputfile(...
    {'*.mat','Files (*.mat)'},'Select a file');
delete(f_dummy); % delete the dummy figure
figure(app.UIFigure);
% If user cancels
if filename == 0
    return
end
close(progdlg);
Filename = [pathname,filename];

%----------------------------------------------------------------------
% Setup ADMESH input
%----------------------------------------------------------------------
PTS = app.PTS;
PTS.Poly(2:end) = [];

k = 0;
for i = 1 : length(Constraints)

    x = Constraints(i).xy(:,1);
    y = Constraints(i).xy(:,2);

    k = k + 1;
    PTS.Constraints(k).num   = -abs(Constraints(i).num);
    PTS.Constraints(k).xy    = [x,y];
    PTS.Constraints(k).type  = Constraints(i).type;
    PTS.Constraints(k).data  = [];
    PTS.Constraints(k).Kappa = [];
end

if isfield(PTS,'cpplon') && ~isempty(PTS.cpplon)
    PTS = CoordinateConversion(app,PTS,'reverse');
end
save(Filename,'PTS');

app.OutputSaved = 1;
app.SavedFile = Filename;