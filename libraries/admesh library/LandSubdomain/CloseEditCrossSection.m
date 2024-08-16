function CloseEditCrossSection(app)

fig = findobj('Type','Figure','Name','Cross-section');
delete(fig);

h = findobj(app.UIAxes.Children,'tag','edit-cs');
delete(h);