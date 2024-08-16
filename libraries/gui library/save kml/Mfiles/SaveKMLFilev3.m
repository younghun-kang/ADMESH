function SaveKMLFilev3(varargin)
% SaveKMLFile - GUI Callback that saves mesh output to kml file
%
% Syntax:  SaveKMLFile(varargin)
%
% Inputs:
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

%--------------------------- BEGIN CODE -----------------------------------

%--------------------------------------------------------------------------
% Get GUI data
%--------------------------------------------------------------------------

app = varargin{1};

MESH = app.MESH; clear guiH

%--------------------------------------------------------------------------
% Check for variables
%--------------------------------------------------------------------------
if isempty(MESH) % User has not run ADmesh yet
    msg = 'No mesh to save....';
    uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');
    return
end

%--------------------------------------------------------------------------
% Prepare mesh for kml write out. 
%--------------------------------------------------------------------------

% Convert coordinates
if isfield(MESH,'cpplon')
    MESH = Meters2Geo(MESH);
else
    errordlg(['In order to write a kml file the coordinates ',...
        'of the domain needed to be geographic originally.'])
    return
end

% Ask user for file name
% gui.sb.setText('Save file as...')
[file,path] = uiputfile('*.kml','Save As');

% If user cancels
if isequal(file,0) || isequal(path,0)
    % gui.sb.setText('Ready')
    return
end

%--------------------------------------------------------------------------
% Break down each element into line segments
%--------------------------------------------------------------------------
% Assign nodes and connectivity list
p = MESH.Points;
t = MESH.ConnectivityList;

% Get triangulation representation
trep = triangulation(t,p);

% Get edges, sort left to right
E = sort(edges(trep),2);

% Initialize cell for holding boundary conditions
BC = cell(numel(MESH.Constraints),1);

% Remove boundary conditions from edges matrix
for k = 1:numel(MESH.Constraints)
    
    % Need to handle external barriers differently
    if any(MESH.Constraints(k).num == [4 5 24 25])
        
        % 1st Node string
        nStr1 = MESH.Constraints(k).nodeStr(:,1);
        % 2nd Node string
        nStr2 = MESH.Constraints(k).nodeStr(:,2);
        
        % Compete edge list
        ns = [[nStr1(1:end-1),nStr1(2:end)];[nStr2(1:end-1),nStr2(2:end)]];
        
        % Append edge type
        type = ones(size(ns,1),1)*MESH.Constraints(k).num;
        
        % Constraint (append type)
        BC{k} = [sort(ns,2) type];

    else
        
        % Node string
        ns = MESH.Constraints(k).nodeStr;
        
        % Append edge type
        type = ones(size(ns,1)-1,1)*MESH.Constraints(k).num;
        
        % Constraint (append type)
        BC{k} = [sort([ns(1:end-1), ns(2:end)],2) type];
        
    end
    
    E = setdiff(E,BC{k}(:,[1 2]),'rows'); % Remove common points
    
end

% For regular element edges, append nan
E = [E nan(size(E,1),1)];

% Concatenate elements edges
E = [cell2mat(BC); E];

%--------------------------------------------------------------------------
% Determine if file should be broken up
%--------------------------------------------------------------------------
% var = whos('E');
% 
% if var.bytes > 1000 % Break up
%     
%     
%     
%     
%     
% end


%--------------------------------------------------------------------------
% Prepare to write
%--------------------------------------------------------------------------

% Initialize kml header
Header = {...
    '<?xml version="1.0" encoding="UTF-8"?>',...
    '<kml xmlns="http://www.opengis.net/kml/2.2">',...
    '   <Document>',...
    '   <Style id="yellowStyle">',...
    '       <LineStyle>',...
    '           <color>ff00ffff</color>',...
    '           <width>1</width>',...
    '       </LineStyle>',...
    '   </Style>',...
    '   <Style id="redStyle">',...
    '       <LineStyle>',...
    '           <color>ff0000ff</color>',...
    '           <width>2</width>',...
    '       </LineStyle>',...
    '   </Style>',...
    '   <Style id="blueStyle">',...
    '       <LineStyle>',...
    '           <color>ffff0000</color>',...
    '           <width>2</width>',...
    '       </LineStyle>',...
    '   </Style>',... 
    '   <Style id="greenStyle">',...
    '       <LineStyle>',...
    '           <color>ff00ff00</color>',...
    '           <width>2</width>',...
    '       </LineStyle>',...
    '   </Style>',... 
    '   <Style id="blackStyle">',...
    '       <LineStyle>',...
    '           <color>ff000000</color>',...
    '           <width>2</width>',...
    '       </LineStyle>',...
    '   </Style>',...
    '   <Style id="brownStyle">',...
    '       <LineStyle>',...
    '           <color>ff336699</color>',...
    '           <width>2</width>',...
    '       </LineStyle>',...
    '   </Style>'}';

% Initialize beginning & ending block used for each path
bBlock = {...
    '   <Placemark>',...
    '   <name>Element</name>',...
    '   <styleUrl>#yellowStyle</styleUrl>',...
    '   <LineString>',...
    '       <tessellate>1</tessellate>',...
    '       <coordinates>'}';

eBlock = {...
    '       </coordinates>',...
    '   </LineString>',...
    '   </Placemark>'}';

% Initialize kml footer
Footer = {...
    '   </Document>',...
    '</kml>'}';

% Initialze cell for fprintf 
nedges     = size(E,1);
nHeaders   = size(Header,1);
nLines     = 11*nedges;
nFooters   = 2;

textOut = cell(nHeaders + nLines + nFooters,1);

% gui.sb.setText('Preparing mesh for kml write out...')

% Store header
textOut(1:nHeaders) = Header;

% Initialize start and end index for edge write out
s = nHeaders+1;
e = s + 10;

% Determine appropraite frequency for updating status bar
freq = round(nedges/100);

msg = 'Writing kml file...';
progdlg = uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',msg);
progdlg.Value = 1/nedges;

% Loop over each edge and prep for write out
for k = 1:nedges
    
    % Assign (x,y,z) coordinated for edge
    x = [p(E(k,1),1); p(E(k,2),1)];
    y = [p(E(k,1),2); p(E(k,2),2)];
    z = [p(E(k,1),3); p(E(k,2),3)];
    
    % What color are we using?
    if isnan(E(k,3))
        bBlock(3) = {'   <styleUrl>#yellowStyle</styleUrl>'};
    elseif E(k,3) == -1
        bBlock(3) = {'   <styleUrl>#blueStyle</styleUrl>'};
    elseif any(E(k,3) == [0 2 10 12 20 22 30])
        bBlock(3) = {'   <styleUrl>#blackStyle</styleUrl>'}; 
    elseif any(E(k,3) == [1 11 21])
        bBlock(3) = {'   <styleUrl>#greenStyle</styleUrl>'}; 
    elseif any(E(k,3) == [3 13 23 4 5 24 25])
        bBlock(3) = {'   <styleUrl>#brownStyle</styleUrl>'}; 
    elseif E(k,3) == 18
        bBlock(3) = {'   <styleUrl>#redStyle</styleUrl>'};         
    end
    
    textOut(s:e) = [...
        bBlock; ...
        [{sprintf('               %5.10f,%5.10f,%f',x(1),y(1),z(1) )},...
        {sprintf('               %5.10f,%5.10f,%f',x(2),y(2),z(2) )}]';...
        eBlock];
    
    s = e + 1;
    e = s + 10;
    
    if(mod(k,freq ) == 0)
        progdlg.Value = k/nedges;
    end
    
end

% Append footer
textOut(end-1:end) = Footer;

%--------------------------------------------------------------------------
% Open file and write
%--------------------------------------------------------------------------

msg = 'Writing kml file...';
uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');

% Open file
fid = fopen([path file],'w');

fprintf(fid,'%s\n',textOut{:});

% Close file
fclose(fid);

end