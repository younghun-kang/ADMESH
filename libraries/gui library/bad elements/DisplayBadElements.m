function DisplayBadElements(app)
% DisplayBadElements - Toolbar Callback that displays elements with element
% quality below what is desired
%
% Syntax:  DisplayBadElements(guiFig,~,~)
%
% Inputs:
%    guiFig - handle that identifies the figure
%
% Outputs:
%    non
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% August 2013; Last revision: 08-August-2013

%------------- BEGIN CODE --------------


%------------------------------------------------------------------------------
% Get GUI data
%------------------------------------------------------------------------------
pH          = app.UIAxes;
meshPatch   = app.MESH;

%------------------------------------------------------------------------------
% Check if there is anything currently in the plot window
%------------------------------------------------------------------------------
plotItems = get(pH, 'Children');

if isempty(plotItems) || isempty(meshPatch)

    return
    
end

%------------------------------------------------------------------------------
% Ask user for minimum element quality
%------------------------------------------------------------------------------
prompt={[...
    'Enter a minimum element quality value ranging between 0 to 1. 0',...
    ' representing a degenerate triangle and 1 representing an equaliteral',...
    ' triangle. All elements that fall below this threshold will be highlighted',...
    ' in red.']};
name='ADmesh'; numlines=1; defaultanswer={num2str(app.MinEQ)};
options.Resize='off'; options.WindowStyle='normal'; options.Interpreter='tex'; pause(.001)

try
    app.MinEQ = str2double(strtrim(cell2mat(inputdlg(prompt,name,numlines,defaultanswer,options))));
end
drawnow; pause(0.05);  % this innocent line prevents the Matlab hang

% Check range
if app.MinEQ < 0 || app.MinEQ > 1
    
    msg = 'The value you entered was out of the range of [0 , 1]';
    uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');

    return
    
end

%--------------------------------------------------------------------------
% Save gui data
%--------------------------------------------------------------------------
% guidata(fig,gui);

%--------------------------------------------------------------------------
% Display mesh info
%--------------------------------------------------------------------------
DisplayMeshInfo(app)

end