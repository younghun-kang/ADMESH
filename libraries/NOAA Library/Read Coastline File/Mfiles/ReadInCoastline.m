function PTS = ReadInCoastline(guiFig)

%--------------------------------------------------------------------------
% Load GUI Data
%--------------------------------------------------------------------------
guiH = % guidata(guiFig);

pH = guiH.PlotWindow; % Plto axes

%--------------------------------------------------------------------------
% Read in File
%--------------------------------------------------------------------------
[filename, pathname] = uigetfile( '*.txt*', 'Select a file');

% If user cancels
if ~filename; PTS = nan; return; end

set(guiH.SB.ProgressBar, 'Visible','on', 'Indeterminate','on');
statusbar(guiFig, 'Reading in coastline data...'); drawnow

% Open File
fid = fopen([pathname filename],'r');

% Read in (lon,lat) data
Coastline_Data = cell2mat(textscan(fid,'%f %f'));

% Close file
fclose(fid);

%--------------------------------------------------------------------------
% Enter while loop
%--------------------------------------------------------------------------
while 1
    
    %----------------------------------------------------------------------
    % Open dialogue box. Ask user for some info
    %----------------------------------------------------------------------
    statusbar(guiFig, 'Enter spacing tolerance...'); drawnow
    prompt={[...
        'Enter a tolerance or a distance within which the coastline '...
        'segments will be joined:']};
    name='ADmesh'; numlines=1; defaultanswer={'.0001'};
    options.Resize='off'; options.WindowStyle='normal'; 
    options.Interpreter='tex'; pause(.001)
    
    tol = inputdlg(prompt,name,numlines,defaultanswer,options);
    
    drawnow; pause(.005)
    
    % If user cancels
    if isempty(tol); PTS = nan; return; end
    
    % Deal input data
    tol = str2double(strtrim(tol{:})); 
    
    %----------------------------------------------------------------------
    % Seperate mainland and islands based on tolerance
    %----------------------------------------------------------------------
    statusbar(guiFig, 'Seperating mainland and islands based on tolerance...'); drawnow
    
    [Coastline_Data]=join_cst(Coastline_Data,tol);
    
    %----------------------------------------------------------------------
    % Refine coastline and put into edge structure
    %----------------------------------------------------------------------
    statusbar(guiFig, 'Refining coastline and sorting into edge structure...'); drawnow
    
    [PTS] = NOAA_Coastline_Sort(Coastline_Data);
  
    %----------------------------------------------------------------------
    % Plot PTS
    %----------------------------------------------------------------------
    PlotEdgeStructure(PTS,nan,guiFig,.1)
    
    %------------------------------------------------------------------------------------
    % Ask user if they would like to continue
    %------------------------------------------------------------------------------------
    msg = ['Would you like to continue with this '...
        'coastline or adjust the tolerance for a different'...
        ' coastline extraction?'];
    choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'Continue','Re-Try','Quit'},'DefaultOption',1,'Icon','Warning');
    drawnow; pause(.005);  % this innocent line prevents the Matlab hang
    
    switch quest
        case 'Continue'
            
            break
            
        case 'Quit'
            
            PTS = nan;
            break
            
        case 'Re-Try'
            
            % Clear window
            plotItems = get(pH, 'Children');
            
            if ~isempty(plotItems); delete(plotItems); end
            
    end
    
end

end