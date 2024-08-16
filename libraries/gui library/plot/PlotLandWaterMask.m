function PlotLandWaterMask(app)

% if isprop(app,'MainApp')
%     app.UIAxes = app.MainApp.UIAxes;
% end
h = findobj(app.UIAxes,'tag','Land Water Mask');
delete(h);

plot(app.UIAxes,app.pgon,'FaceColor',[0 0.4470 0.7410],'tag','Land Water Mask');

figure(app.MainApp.UIFigure);
figure(app.UIFigure);