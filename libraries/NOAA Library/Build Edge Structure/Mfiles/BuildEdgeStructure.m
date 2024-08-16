function PTS = BuildEdgeStructure(Coastline_Data,pH,gui)

sb = gui.sb;

%--------------------------------------------------------------------------
% Enter while loop
%--------------------------------------------------------------------------
while 1
    
    %----------------------------------------------------------------------
    % Open dialogue box. Ask user for tolerance
    %----------------------------------------------------------------------
    sb.setText('Enter spacing tolerance...')
    
    prompt={[...
        'Enter a tolerance or a distance within which the coastline '...
        'segments will be joined:']};
    name='ADmesh'; numlines=1; defaultanswer={'.0001'};
    options.Resize='off'; options.WindowStyle='normal';
    options.Interpreter='tex'; pause(.001)
    
    tol = inputdlg(prompt,name,numlines,defaultanswer,options);
    
    drawnow; pause(.005)
    
    %----------------------------------------------------------------------
    % Check user input
    %----------------------------------------------------------------------
    
    % If user cancels
    if isempty(tol); PTS = []; sb.setText('Ready'); return; end
    
    % Check is numeric string
    if any(isstrprop(tol{:}, 'alpha'))
        sb.setText('Ready')
        errordlg('The tolerance input is not numeric.','Error');
        PTS = []; return;
    end
    
    % Deal input data
    tol = str2double(strtrim(tol{:})); 
    
    %----------------------------------------------------------------------
    % Seperate mainland and islands based on tolerance
    %----------------------------------------------------------------------
    sb.setText('Seperating mainland and islands based on tolerance...')
    
    [Coastline_Data]=join_cst(Coastline_Data,tol);
    
    %----------------------------------------------------------------------
    % Refine coastline and put into edge structure
    %----------------------------------------------------------------------
    sb.setText('Refining coastline and sorting into edge structure...')

    [PTS] = NOAA_Coastline_Sort(Coastline_Data);
        
    %----------------------------------------------------------------------
    % Convert edge structure to meters
    %----------------------------------------------------------------------
    [PTS,PTS.cpplon,PTS.cpplat] = Geo2Meters(PTS);
    
    %----------------------------------------------------------------------
    % Plot PTS
    %----------------------------------------------------------------------
    PlotEdgeStructure(PTS,gui,sb,pH,.1)
    
    %---------------------------------------------------------------------------
    % Ask user if they would like to continue
    %---------------------------------------------------------------------------
    sb.setText('Continue?')
    
    msg = ['Would you like to continue with this '...
        'coastline or adjust the tolerance for a different'...
        ' coastline extraction?'];
    choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'Continue','Re-Try','Quit'},'DefaultOption',1,'Icon','Warning');
    drawnow; pause(.005);  % this innocent line prevents the Matlab hang
    
    switch quest
        case 'Continue'
            
            break % Break while loop and continue
            
        case 'Quit'
            
            PTS = [];
            break
            
        case 'Re-Try'
            
            % Clear window
            plotItems = get(pH, 'Children');
            
            if ~isempty(plotItems); delete(plotItems); end
            
    end
    
end

end