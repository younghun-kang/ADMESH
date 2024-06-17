function MinAreaChangedExtractChannelsWaterApp(app)

if isempty(app.pgon)
    return;
end

msg = 'You may want to re-create Land-Water mask after changing "Min Area".';
choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
    'Options',{'Re-apply','Ignore'},'DefaultOption',1,'Icon','warning');

if strcmp(choice,'Ignore')
    return;
end

CreateLandWaterMask(app);

