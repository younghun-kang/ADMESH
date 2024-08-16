function [h_gradlimit,H] = GradientLimiting(x,h,g)

%--------------------------------------------------------------------------
% Set constant parameters
%--------------------------------------------------------------------------
MaxIter = 1e3;
IterTol = 1e-6;
dt = min(diff(x))/2;

h_old = h(x);

id_stand = 2:length(x)-1;
inv_dxf = 1 ./ (x(id_stand+1) - x(id_stand));
inv_dxb = 1 ./ (x(id_stand) - x(id_stand-1));

nIter = 0;
H(nIter+1,:) = h_old;
while 1
    %----------------------------------------------------------------------
    % Compute gradient using Persson's formula
    %----------------------------------------------------------------------
    fd = (h_old(id_stand+1) - h_old(id_stand)) .* inv_dxf;
    bd = (h_old(id_stand) - h_old(id_stand-1)) .* inv_dxb;
    
    FD = [(h_old(2) - h_old(1)) ./ (x(2) - x(1)); fd(:); 0];
    BD = [0; bd(:); (h_old(end) - h_old(end-1)) ./ (x(end) - x(end-1))];
    
    gradh1 = sqrt( max(BD,0).^2 + min(FD,0).^2 );
    
    %----------------------------------------------------------------------
    % Compute RHS and update mesh size
    %----------------------------------------------------------------------
    RHS = (min(gradh1,g) - gradh1);
    h_new = h_old + dt*RHS;
    h_diff = h_new - h_old;
    h_old = h_new;
    
    %----------------------------------------------------------------------
    % Check if steady state
    %----------------------------------------------------------------------
    nIter = nIter + 1;
    H(nIter+1,:) = h_old;
    if max(gradh1) - g < IterTol
        break;
    end
    
    if nIter > MaxIter
        warning([...
            'Potential issues are in GradientLimiting program.',...
            'Program ended as maximum iteration reached.']);
        break;
    end
    
end

h_gradlimit = @(y) interp1(x,h_old,y);


end



