function [output,status] = CoordinateConversion(app,input,mode,varargin)

output = input;
status = 0;

if isempty(input)    
    return;
end

if length(varargin) == 2
    cpplon = varargin{1};
    cpplat = varargin{2};
elseif isempty(varargin)
    try
        cpplon = input.cpplon;
        cpplat = input.cpplat;
    catch
        cpplon = [];
        cpplat = [];
    end
end

progdlg = [];

switch lower(mode)
    case 'forward' % Convert to XY (meters)
        progdlg = uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',...
            'Converting to Cartersian coordinates...','Indeterminate','on');
        if ~isempty(cpplon)
            output = Geo2Cart(input,cpplon,cpplat);
            status = 1;
        else
            [output,output.cpplon,output.cpplat] = Geo2Cart(input);
            status = 1;
        end

    case 'reverse' % Convert to lat/lon
        progdlg = uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',...
            'Converting to Geographic coordinates...','Indeterminate','on');
        output = Cart2Geo(input);
        status = 1;

    case 'auto'
        if isfield(input,'cpplon') && ~isempty(input.cpplon)
            if isfield(input,'Poly') % PTS structure input
                progdlg = uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',...
                    'Converting to Cartersian coordinates...','Indeterminate','on');
                output = Geo2Cart(input,input.cpplon,input.cpplat);
                status = 1;
            elseif isfield(input,'Points') % MESH
                progdlg = uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',...
                    'Converting to Geographic coordinates...','Indeterminate','on');
                output = Cart2Geo(input,input.cpplon,input.cpplat);
                status = 1;
            end
        else
            if isfield(input,'Poly') % PTS structure input
                x = vertcat(input.Poly.x);
                y = vertcat(input.Poly.y);
                x = x(~isnan(x));
                y = y(~isnan(y));

            elseif isfield(input,'Points') || isprop(input,'Points')% MESH or xyzFun
                x = input.Points(:,1);
                y = input.Points(:,2);

            else
                msg = 'Unknown format of input is passed in CoordinateConversion';
                uiconfirm(app.UIFigure,msg,'ADMESH',...
                    'Options',{'OK'},'DefaultOption',1,'Icon','Error');
            end

            if all(x < 180 & x > -180 & y < 90 & y > -90)

                msg = ['The input may be in geographic coordinate system (in lat/lon). '...
                    'Do you want to convert it into planar coordinate system (in meters)? ',...
                    'If you are unsure, please check the scalebar at the right bottom of the figure.'];
                choice = uiconfirm(app.UIFigure,msg,'Coordinate conversion',...
                    'Options',{'Convert coordinates','Do not convert coordinates'},'DefaultOption',1,'Icon','Warning');

                if strcmpi(choice,'Convert coordinates')
                    % Convert to XY (meters)
                    progdlg = uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',...
                        'Converting to Cartersian coordinates...','Indeterminate','on');
                    [output,output.cpplon,output.cpplat] = Geo2Cart(input);
                    status = 1;
                end

            end
        end

    otherwise
        msg = 'The mode must be one of "forward", "reverse", and "auto".';
        uiconfirm(app.UIFigure,msg,'ADMESH',...
            'Options',{'OK'},'DefaultOption',1,'Icon','Error');

end

close(progdlg);


