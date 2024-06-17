function varargout = Geo2Cart(varargin)
% Geo2Meters - 
%
% Syntax:  varargout = Geo2Meters(varargin)
%
% Other m-files required: READ_14
% Subfunctions: none
% MAT-files required: none
%
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
%---------------------------- BEGIN CODE --------------------------------------

%------------------------------------------------------------------------------
% Convert lat/lon to meters for PTS data structure
%------------------------------------------------------------------------------
if nargin == 1 
    
    if isfield(varargin{1},'Poly') % PTS structure input
        PTS = varargin{1};

        % Radius of earth
        r = 6378206;

        cpplon = (pi/180)*(PTS.Poly(1).x(1)+1.2345678);
        cpplat = (pi/180)*(PTS.Poly(1).y(1)+1.2345678);

        for k = 1:numel(PTS.Poly)

            PTS.Poly(k).x = r.*((pi/180).*PTS.Poly(k).x - cpplon).*cos(cpplat);
            PTS.Poly(k).y = r.*((pi/180).*PTS.Poly(k).y);

        end

        % Convert constraints
        if isfield(PTS,'Constraints')
            for k = 1:numel(PTS.Constraints)

                x = PTS.Constraints(k).xy(:,1);
                y = PTS.Constraints(k).xy(:,2);

                x = r.*((pi/180).*x - cpplon).*cos(cpplat);
                y = r.*((pi/180).*y);

                PTS.Constraints(k).xy = [x y];

            end
        end

        varargout{1} = PTS;
        varargout{2} = cpplon;
        varargout{3} = cpplat;

    elseif isfield(varargin{1},'Points') % MESH

        MESH = varargin{1};

        lon = MESH.Points(:,1);
        lat = MESH.Points(:,2);

        % Radius of earth
        r = 6378206;

        cpplon = (pi/180)*(lon(1)+1.2345678);
        cpplat = (pi/180)*(lat(1)+1.2345678);

        x = r.*((pi/180).*lon - cpplon).*cos(cpplat);
        y = r.*((pi/180).*lat);

        MESH.Points(:,[1 2]) = [x,y];

        varargout{1} = MESH;
        varargout{2} = cpplon;
        varargout{3} = cpplat;

    elseif isa(varargin{1},'scatteredInterpolant') % xyzFun (scatteredInterpolant)

        xyzFun = varargin{1};

        lon = xyzFun.Points(:,1);
        lat = xyzFun.Points(:,2);

        % Radius of earth
        r = 6378206;

        cpplon = (pi/180)*(lon(1)+1.2345678);
        cpplat = (pi/180)*(lat(1)+1.2345678);

        x = r.*((pi/180).*lon - cpplon).*cos(cpplat);
        y = r.*((pi/180).*lat);

        xyzFun.Points(:,[1 2]) = [x,y];

        varargout{1} = xyzFun;

    elseif isa(varargin{1},'griddedInterpolant') % xyzFun (griddedInterpolant)

        xyzFun = varargin{1};

        lon = xyzFun.GridVectors{1};
        lat = xyzFun.GridVectors{2};

        % Radius of earth
        r = 6378206;

        cpplon = (pi/180)*(lon(1)+1.2345678);
        cpplat = (pi/180)*(lat(1)+1.2345678);

        x = r.*((pi/180).*lon - cpplon).*cos(cpplat);
        y = r.*((pi/180).*lat);

        xyzFun.GridVectors{1} = x;
        xyzFun.GridVectors{2} = y;

        varargout{1} = xyzFun;

    end
    
%------------------------------------------------------------------------------
% Convert lat/lon to meters for lon/lat input
%------------------------------------------------------------------------------
elseif nargin == 2 
    
    lon = varargin{1};
    lat = varargin{2};
    
    % Radius of earth
    r = 6378206;
    
    cpplon = (pi/180)*(lon(1)+1.2345678);
    cpplat = (pi/180)*(lat(1)+1.2345678);
    
    
    x = r.*((pi/180).*lon - cpplon).*cos(cpplat);
    y = r.*((pi/180).*lat);
    
    
    varargout{1} = x;
    varargout{2} = y;
    varargout{3} = cpplon;
    varargout{4} = cpplat;    
        
%------------------------------------------------------------------------------
% Convert lat/lon to meters for PTS data structure with cpplon/cpplat
%------------------------------------------------------------------------------
elseif nargin == 3
    
    cpplon = varargin{2};
    cpplat = varargin{3};

    if isfield(varargin{1},'Poly') % PTS structure input
        PTS = varargin{1};

        % Radius of earth
        r = 6378206;

        for k = 1:numel(PTS.Poly)

            PTS.Poly(k).x = r.*((pi/180).*PTS.Poly(k).x - cpplon).*cos(cpplat);
            PTS.Poly(k).y = r.*((pi/180).*PTS.Poly(k).y);

        end

        % Convert constraints
        if isfield(PTS,'Constraints')
            for k = 1:numel(PTS.Constraints)

                x = PTS.Constraints(k).xy(:,1);
                y = PTS.Constraints(k).xy(:,2);

                x = r.*((pi/180).*x - cpplon).*cos(cpplat);
                y = r.*((pi/180).*y);

                PTS.Constraints(k).xy = [x y];

            end
        end

        varargout{1} = PTS;

    elseif isfield(varargin{1},'Points') % MESH

        MESH = varargin{1};

        lon = MESH.Points(:,1);
        lat = MESH.Points(:,2);

        % Radius of earth
        r = 6378206;

        x = r.*((pi/180).*lon - cpplon).*cos(cpplat);
        y = r.*((pi/180).*lat);

        MESH.Points(:,[1 2]) = [x,y];

        varargout{1} = MESH;

    elseif isa(varargin{1},'scatteredInterpolant') % xyzFun (scatteredInterpolant)

        xyzFun = varargin{1};

        lon = xyzFun.Points(:,1);
        lat = xyzFun.Points(:,2);

        % Radius of earth
        r = 6378206;

        x = r.*((pi/180).*lon - cpplon).*cos(cpplat);
        y = r.*((pi/180).*lat);

        xyzFun.Points(:,[1 2]) = [x,y];

        varargout{1} = xyzFun;

    elseif isa(varargin{1},'griddedInterpolant') % xyzFun (griddedInterpolant)

        xyzFun = varargin{1};

        lon = xyzFun.GridVectors{1};
        lat = xyzFun.GridVectors{2};

        % Radius of earth
        r = 6378206;

        x = r.*((pi/180).*lon - cpplon).*cos(cpplat);
        y = r.*((pi/180).*lat);

        xyzFun.GridVectors{1} = x;
        xyzFun.GridVectors{2} = y;

        varargout{1} = xyzFun;

    end

%------------------------------------------------------------------------------
% Convert lat/lon to meters for bathymetry with cpplon/cpplat
%------------------------------------------------------------------------------    
elseif nargin == 4
    
    lon     = varargin{1};
    lat     = varargin{2};
    cpplon  = varargin{3};
    cpplat  = varargin{4};
    
    % Radius of earth
    r = 6378206;
    
    x = r.*((pi/180).*lon - cpplon).*cos(cpplat);
    y = r.*(pi/180).*lat;
    
    varargout{1} = x;
    varargout{2} = y;
    
end


end