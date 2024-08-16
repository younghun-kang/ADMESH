function CloseADMESH(app)

msg = 'Are you sure you want to close ADMESH?';
choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
    'Options',{'Yes','No'},'DefaultOption',2,'Icon','Warning');
drawnow; pause(.005);

switch choice
    
    case 'Yes'
        % Close subapps
        for i = 1 : length(app.subapps)
            subappnames = fieldnames(app.subapps);
            subapp = app.subapps.(subappnames{i});
            if isvalid(subapp)
                delete(subapp);
            end
        end
        % Close ADMESH app
        delete(app);
        % Remove libraries from path
        warning_status = warning;
        warning('off');
        rmpath(genpath('libraries'));
        warning(warning_status);
    case 'No'
        
        
end



end