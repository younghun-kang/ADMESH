function RetrieveChannelBanksOnConstraints(app)

PTS = app.PTS;
option = app.RetrieveBankEditField.Value;
Mesh1D = [];
k = 0;
for i = 1 : length(PTS.Constraints)
    if PTS.Constraints(i).num == 18
        k = k + 1;
        Mesh1D(k).X = PTS.Constraints(i).xy(:,1);
        Mesh1D(k).Y = PTS.Constraints(i).xy(:,2);
    end
end
%%
m2deg = km2deg(1e-3);
xyzFun = app.xyzFun;
CoordinateSystem = app.CoordinateSystemDropDown.Value;
if strcmpi(CoordinateSystem,'Projected (m)')
    MaxChannelWidth = app.MaxChannelWidthEditField.Value;
    delta = 1;
    unit = 'meter';
elseif strcmpi(CoordinateSystem,'Unprojected (decimal degree)')
    MaxChannelWidth = app.MaxChannelWidthEditField.Value*m2deg;
    delta = 1*m2deg;
    unit = 'degree';
end


CrossSection = app.CrossSection;

app.ProgressBarButton.Text = 'Retrieving banks...'; drawnow;
for i = 1 : length(Mesh1D)

    Centerline = [Mesh1D(i).X,Mesh1D(i).Y];

    PARAMS_bank = RetrieveBankline(xyzFun,Centerline,MaxChannelWidth,delta,option,unit);
    app.ProgressBarButton.Text = sprintf('Retrieving banks... (%d/%d)',i,length(Mesh1D));
    drawnow;

    if isstruct(CrossSection{i})
        CrossSection{i} = struct2table(CrossSection{i});
    end
    
    CrossSection{i}.db = PARAMS_bank.d;
    CrossSection{i}.zb = PARAMS_bank.z;
    CrossSection{i}.cc = PARAMS_bank.c;
    % CrossSection{i} = table2struct(CrossSection{i});
end

app.CrossSection = CrossSection;

PlotBanklines(app);

app.ProgressBarButton.Text = 'Ready'; drawnow;



