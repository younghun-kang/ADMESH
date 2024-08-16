function PlotConstraints(app,pH)
% PlotConstraints - Plots constraints defined in PTS
%
% Syntax:  PlotConstraints(PTS,per)
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% August 2013; Last revision: 12-August-2013

%--------------------------- BEGIN CODE ---------------------------------------

%------------------------------------------------------------------------------
% Initialize linewidth
%------------------------------------------------------------------------------
LineWidth = 1.5;

%------------------------------------------------------------------------------
% Tell ADMESH where to plot
%------------------------------------------------------------------------------
% axes(pH); 
hold(pH,'on');

%------------------------------------------------------------------------------
% Define a context menu; it is not attached to anything
%------------------------------------------------------------------------------
% hcmenu = uicontextmenu;
% 
% % Define the context menu items and install their callbacks for constraints
% uimenu(hcmenu,'Label','Remove Boundary Constraint','Callback',@RemoveBC);

%------------------------------------------------------------------------------
% Plot constraints
%------------------------------------------------------------------------------
PTS = app.PTS;
if isfield(PTS,'Constraints') && ~isempty(PTS.Constraints)
    Constraints = struct2table(PTS.Constraints,'AsArray',true);

    % Internal Constraints
    PlotConstraintsSub([4 24 5 25],'r',LineWidth,'Internal Constraint');

    % External Constraints
    PlotConstraintsSub([3 13 23],'r',LineWidth,'External Constraint');

    % Open Ocean
    PlotConstraintsSub(-1,'b',LineWidth,'Open Ocean');

    % Line
    PlotConstraintsSub([18 17 19],'r',LineWidth,'Line Constraint');
    PlotConstraintsSub(-[18 17 19],'r',LineWidth,'Line Constraint - dummy');

end

    function PlotConstraintsSub(NumConstraint,Color,LineWidth,Tag)
        ix = find(ismember([PTS.Constraints.num],NumConstraint));

        if isempty(ix)
            return;
        end

        p = Constraints.xy(ix);
        p = cellfun(@(x) [x; nan(1,size(x,2))],p,'UniformOutput',0);
        p = vertcat(p{:});

        h = line(p(:,1),p(:,2),'Parent',pH);

        set(h,...
            'Color',Color,...
            'LineWidth',LineWidth,...
            'tag',Tag,...
            'UserData',ix)

        if contains(Tag,'- dummy')
            set(h,'LineStyle','--');
        end

        if strcmpi(Tag,'Line Constraint')
            h.ButtonDownFcn = @(hObj,event)EditCrossSection(app,hObj,event);
        end

    end

end
