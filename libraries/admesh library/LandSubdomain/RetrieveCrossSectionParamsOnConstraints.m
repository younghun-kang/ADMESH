function RetrieveCrossSectionParamsOnConstraints(app)

CS = app.CrossSection;

app.ProgressBarButton.Text = 'Retrieving cross section parameters...'; drawnow;
for i = 1 : length(CS)
    CSParams = RetrieveCrossSectionParameters(app,i);

    CS{i}.TW = CSParams.TW;
    CS{i}.BW = CSParams.BW;
    CS{i}.H = CSParams.H;
    CS{i}.CSC = CSParams.CSC;
    CS{i}.ERR = CSParams.ERR;
    CS{i}.BEL = CSParams.BEL;
    
    app.ProgressBarButton.Text = sprintf('Retrieving cross section parameters... (%d/%d)',i,length(CS));
end

app.CrossSection = CS;

app.ProgressBarButton.Text = 'Ready'; drawnow;



