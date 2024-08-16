function InterpolateCrossSectionsOnConstraints(app)

PTS = app.PTS;

Channels = cell(length(PTS.Constraints),1);
k = 0;
for i = 1 : length(PTS.Constraints)
    if PTS.Constraints(i).num == 18
        k = k + 1;
        Channels{k} = PTS.Constraints(i).xy(:,[1 2]);
    end
end
Channels(k+1:end) = [];

%%
xyzFun = app.xyzFun;
CS = cell(length(Channels),1);

app.ProgressBarButton.Text = 'Retrieving cross sections...'; drawnow;

for i = 1 : length(Channels)
    Centerline = Channels{i};
    
    CS{i} = InterpolateCrossSections(xyzFun,Centerline);
    
    app.ProgressBarButton.Text = sprintf('Retrieving cross sections... (%d/%d)',i,length(Channels));
    drawnow;
end

app.CrossSection = CS;

app.ProgressBarButton.Text = 'Ready'; drawnow;



