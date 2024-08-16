function h0 = CurvatureFunction(h0,D,gradD,X,Y,K,g,hmax,hmin,Settings,UIFigure,PTS)
% curvature - Computes initial mesh size h_curve based on boundary
% curvature
%
% Syntax:  h_curve = curvature(D,gradD,X,Y,K,g,delta,guiFig)
%
% Inputs:
%    D - Distance function computed by the distance function.
%    X - x-coordinates to rectangular grid
%    Y - y-coordinates to rectangular grid
%    delta - Grid spacing
%    g - grading requirement
%    K - number of elements per radian
%    guiFig - handle that identifies the figure
%
% Outputs:
%    h_curve - mesh size based on curvature
%
% Other m-files required: none
% Subfunctions: meshgrid
% MAT-files required: none

% Author: Dustin West
% Equation obtained from Mesh Generation for Implicit Geometries by Per-Olof Persson
% The Ohio State University
% email address: dww.425@gmail.com
% August 2013; Last revision: 10-August-2013

%------------- BEGIN CODE --------------

%------------------------------------------------------------------------------
% Curvature Status
%------------------------------------------------------------------------------
Status = Settings.K.Status;

%------------------------------------------------------------------------------
% Compute Curavture, If On.
%------------------------------------------------------------------------------
if strcmp(Status,'On') 
    
    msg = 'Computing Boundary Curvature...';
    uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');

    % Compute the magnitude of the gradient
    m = ( sqrt(gradD.x.^2 + gradD.y.^2) );
    
    % Compute the curvature using the divergence formula ( See Ron
    % Goldman,Per-Olof Persson)
    kappa = abs(divergence(X,Y,gradD.x./m,gradD.y./m));
    
    % Find points that are within abs(D) <= 2*delta
    I = abs(D) <= 2*hmin;
    
    % Initialize the curvature mesh size function
    h_curve = ones(size(D)).*hmax;
    
    % Specify the number of elements per radian
    %K = K/pi;
    
    % Compute the mesh size for the narrow band around the boundary. This
    % equation factors in the distance function and grading to best obatin
    % actual curvature values since we are using the background mesh nodes,
    h_curve(I) = abs( (1 + kappa(I).*abs(D(I)))./ (K.*kappa(I)) ) - g.*D(I);
    
    
    %......................................................................
    % 2021-04-07 Younghun: not sure why below blocks needed... commented out
    %......................................................................
    %...Younghun added (compute target size from 1D curvature)
    if ~isempty(PTS.Constraints) && isfield(PTS.Constraints,'Kappa')
    p_const = vertcat(PTS.Constraints(:).xy);
    FLkappa = vertcat(PTS.Constraints(:).Kappa);
    temp = struct2table(PTS.Constraints,'AsArray',1);
    temp = temp.xy;
    [~,~,D1] = VectorDistanceTransform(temp,[X(I),Y(I)]);
    [id_const,Dist2const] = knnsearch(p_const,[X(I),Y(I)],'k',2);
    
    I2 = find(I);
%     id_XYI = Dist2const(:,1) <= hmax;
    id_XYI = D1 <= 2*hmin;
%     id_XYI = true(size(I2));
    w = 1./Dist2const(id_XYI,:);
    kappa(I2(id_XYI)) = sum(FLkappa(id_const(id_XYI,:)).*w,2)./sum(w,2);
    h_curve(I2(id_XYI)) = abs(1./(K.*kappa(I2(id_XYI)))) - 0*g.*D(I2(id_XYI));
    end
    %......................................................................
    % Younghun added: Enforce target mesh size as 1D element size
    %......................................................................
    if ~isempty(PTS.Constraints) && Settings.DummyConstraint == 0
    p_const = vertcat(PTS.Constraints(:).xy);
    mp_const = (p_const(1:end-1,:) + p_const(2:end,:))/2;
    dp_const = sqrt( (p_const(1:end-1,1) - p_const(2:end,1)).^2 +...
                  (p_const(1:end-1,2) - p_const(2:end,2)).^2);
              
    temp = struct2table(PTS.Constraints,'AsArray',1);
    temp = temp.xy;
    [~,~,D1] = VectorDistanceTransform(temp,mp_const);
    mp_const = mp_const(D1 < 2*hmin,:);
    dp_const = dp_const(D1 < 2*hmin);
    
    [id_const,Dist2const] = knnsearch(mp_const,[X(I),Y(I)],'k',2);

    [~,~,D1] = VectorDistanceTransform(temp,[X(I),Y(I)]);
%     id_XYI = Dist2const(:,1) <= hmax;
    I2 = find(I);
    id_XYI = D1 <= 2*hmin;

    w = 1./Dist2const(id_XYI,:);
    h_curve(I2(id_XYI)) = sum(dp_const(id_const(id_XYI,:)).*w,2)./sum(w,2);
    end
    
    % Enforce boundary conditions
    h_curve(h_curve < hmin) = hmin;
    h_curve(h_curve > hmax) = hmax;
    
    % Compare initial conditions and save
    h0 = min(h_curve, h0);
        
end

% FLxy = vertcat(PTS.Constraints(:).xy);
% FLkappa = vertcat(PTS.Constraints(:).Kappa);
% [iFL2XY,dXY2FL] = knnsearch(FLxy,[X(I),Y(I)],'k',1);
% iXY2FL = dXY2FL <= 2*hmin;
% I2 = find(I);
% I2 = I2(iXY2FL);
% figure; plot(X(I2),Y(I2),'.');
% figure; myScatter3(X(I2),Y(I2),FLkappa(iFL2XY(iXY2FL)),10);
% figure; myScatter3(FLxy(iFL2XY,1),FLxy(iFL2XY,2),FLkappa(iFL2XY),10);
% 
% FLxy = vertcat(PTS.Constraints(:).xy);
% FLkappa = vertcat(PTS.Constraints(:).Kappa);
% [iXY2FL,dXY2FL] = knnsearch([X(I),Y(I)],FLxy);
% I2 = find(I);
% I2 = I2(iXY2FL);
% figure; plot(X(I2),Y(I2),'.');
% figure; myScatter3(X(I2),Y(I2),FLkappa,10);

end
















