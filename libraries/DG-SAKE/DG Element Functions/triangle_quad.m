function [x,y,w] = triangle_quad(p,varargin)

%--------------------------------------------------------------------------
%
%  [X,Y,W] = triangle_quad(p)
%  [X,Y,W] = triangle_quad(p,T)
%  [X,Y,W] = triangle_quad(...,Type)
%
%  Returns the weights and points of a p-degree numerical integration form-
%  ula for the triangle.  
%
%--------------------------------------------------------------------------
%
%  Input:
%  ------
%         p:  degree of the polynomial (required).
%         T:  triangle over which the integration is performed given in the 
%             form [x1,y1;x2,y2;x3,y3], where x_i, y_i is the i-th vertice 
%             of the triangle (optional). The default is the reference 
%             triangle defined by [-1,-1; 1,-1; -1,1].
%      Type:  type of numerical integration formula, 'product' or 
%             'nonproduct' (optional). The default is to return the formula 
%             with the fewest points.
%
%  Output:
%  -------
%    X,Y:  locations of the quadrature points
%    W:    corresponding weights
%
%--------------------------------------------------------------------------
%
%  NOTES: (1) A p-degree numerical integration formula integrates all 
%             complete polynomials of degree <= p exactly. 
%         (2) Only so-called PI formulas are included here, i.e., formulas 
%             with all weights 'P'ositive and all points 'I'nside the
%             triangle.
%         (2) Attempts have been made to use the most efficient PI formulas 
%             available, i.e., those that require the fewest points for a 
%             given degree. Relevant works are noted for each degree with 
%             the complete reference listed at the bottom. 
%
%--------------------------------------------------------------------------
%
%  Written by Ethan Kubatko, October 26, 2011
%
%--------------------------------------------------------------------------

% Validate data passed to function
%----------------------------------

vp = @(x)validateattributes(x,{'numeric'},{'scalar','integer','>=',0});
vT = @(x)validateattributes(x,{'numeric'},{'size',[3,2]});
vq = @(x)any(strcmpi(x,{'product', 'nonproduct'}));
ip = inputParser;
ip.addRequired('p',vp); 
if (nargin == 2 && (strcmp(class(varargin{1}),'double') == 1))|nargin == 3
    ip.addOptional('T',[ -1,-1; 1,-1; -1, 1 ],vT);  
    ip.addOptional('Type','nonproduct',vq);     
else
    ip.addOptional('Type','nonproduct',vq);    
    ip.addOptional('T',[ -1,-1; 1,-1; -1, 1 ],vT);
end
ip.parse(p,varargin{:});
ip.Results; T = ip.Results.T; Type = ip.Results.Type;

% Assign vertices of triangle and compute area
%----------------------------------------------

a = [T(1,1),T(1,2)]; b = [T(2,1),T(2,2)]; c = [T(3,1),T(3,2)];
area = 1/2*abs(det([a(1), a(2), 1; b(1), b(2), 1; c(1), c(2), 1]));

if strcmpi('nonproduct',Type) == 1 & p < 11

%**************************************************************************
% p = 1                                                (Lower bound = 1 pt)            
%**************************************************************************
    if p == 0 || p == 1    
    
% 1 Point, PI, Optimal (Stroud, verified)
%--------------------------------------------------------------------------
        m = [ 1 ]; b1 = [ 1/3 ]; b2 = [ 1/3 ]; w1 = [ 1 ];
    
%**************************************************************************
% p = 2                                               (Lower bound = 3 pts)
%**************************************************************************   
    elseif p == 2

% 3 Points, PI, Optimal (Stroud, verified)
%--------------------------------------------------------------------------
        %m = [ 3 ]; b1 = [ 1/6 ]; b2 = [ 1/6 ]; w1 = [ 1/3 ];
        m = [ 3 ]; b1 = [  0  ]; b2 = [ 1/2 ]; w1 = [ 1/3 ];
    
%**************************************************************************
% p = 3                                               (Lower bound = 4 pts)
%**************************************************************************
    elseif p == 3    
    
% 4 Points, PI, Optimal (Hillion, verfied)
%--------------------------------------------------------------------------    
%         m  = [ 1, 1, 1, 1]; 
%         b1 = [   1/3,   1/5,   3/5,   1/5 ];
%         b2 = [   1/3,   3/5,   1/5,   1/5 ];
%         w1 = [ -9/16, 25/48, 25/48, 25/48 ];
          m  = [ 1, 1, 1, 1 ];
          b1 = [ 0.666390246, 0.178558728, 0.280019915, 0.075031109 ];
          b2 = [ 0.178558728, 0.666390246, 0.075031109, 0.280019915 ];  
          w1 = [ 0.159020691, 0.159020691, 0.090979309, 0.090979309 ]*2;
    
%**************************************************************************
% p = 4                                               (Lower bound = 6 pts)
%**************************************************************************  
    elseif p == 4
    
% 6 Points, PI, Optimal (Dunavant, verfied)
%--------------------------------------------------------------------------    
        m  = [ 3, 3 ];
        b1 = [ 0.108103018168070, 0.816847572980459 ];
        b2 = [ 0.445948490915965, 0.091576213509771 ];
        w1 = [ 0.223381589678011, 0.109951743655322 ];    
    
%**************************************************************************
% p = 5                                               (Lower bound = 7 pts)
%**************************************************************************   
    elseif p == 5
    
% 7 Points, PI, Optimal (Lyness, verified)
%--------------------------------------------------------------------------      
        m  = [ 1, 3, 3 ];
        b1 = [  1/3, (  9 + 2*sqrt(15))/21,   (9   - 2*sqrt(15))/21   ];
        b2 = [  1/3, (  6 -   sqrt(15))/21,   (6   +   sqrt(15))/21   ];
        w1 = [ 9/40, (155 -   sqrt(15))/1200, (155 +   sqrt(15))/1200 ];
    
%**************************************************************************
% p = 6                                              (Lower bound = 10 pts)
%**************************************************************************     
    elseif p == 6
    
% 12 Points, PI (Dunavant, verified)
%--------------------------------------------------------------------------    
        m  = [ 3, 3, 6 ];
        b1 = [ 0.501426509658179, 0.873821971016996, 0.053145049844817 ]; 
        b2 = [ 0.249286745170910, 0.063089014491502, 0.310352451033784 ];
        w1 = [ 0.116786275726379, 0.050844906370207, 0.082851075618374 ];        
        
% 10 Points, PO, Optimal (Griener)
%--------------------------------------------------------------------------         
% Griener Example 2 (1 point outside)
%         m  = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ];
%         b1 = [  0.080569513867192,  0.080569513867192, 0.075915879482650, ...
%                 0.075915879482650,  0.384332364623422, 0.367209760982061, ...
%                 0.367209760982061,  1.345107271456660, 0.753927806173285, ...
%                 0.753927806173285 ];
%         b2 = [  0.296812064432131,  0.622618421700677, 0.058482749470777, ...
%                 0.865601371046573,  0.307833817688289, 0.068874709963361, ...
%                 0.563915529054578, -0.172553635728330, 0.048128941614260, ...
%                 0.197943252212456 ];        
%         w1 = [  0.249865692758219,  0.249865692758219, 0.114737978272426, ...
%                 0.114737978272426,  0.420752307674919, 0.250959273633634, ...
%                 0.250959273633634,  0.000855673313626, 0.173633064841449, ...
%                 0.173633064841449 ]; 
% Griener Example 1 (2 points outside)
%         m  = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ];
%         b1 = [ -0.279887910429196, -0.279887910429196, 0.147081984377125, ...
%                 0.147081984377125,  0.521445509916158, 0.521445509916158, ...
%                 0.055642733460735,  0.055642733460735, 0.829333677659378, ...
%                 0.292618906660071 ];
%         b2 = [ -0.008791997144206,  1.288679907573402, 0.064688806056011, ...
%                 0.788229209566864,  0.076825716849082, 0.401728773234760, ...
%                 0.295025119369415,  0.649332147169850, 0.085333161170311, ...
%                 0.353690546669964 ]; 
%         w1 = [  0.001030857518569,  0.001030857518569, 0.209590626265694, ...
%                 0.209590626265694,  0.297102218602087, 0.297102218602087, ...
%                 0.189954478355135,  0.189954478355135, 0.181987047180899, ...
%                 0.422656591336133 ];
%**************************************************************************
% p = 7                                              (Lower bound = 12 pts)
%**************************************************************************
    elseif p == 7

% 12 Points, PI, Optimal (Gatermann, verified) 
%-------------------------------------------------------------------------- 
        m  = [ 3, 3, 3, 3];
        b1 = [ 0.6238226509439084e-1, 0.5522545665692000e-1, ...
               0.3432430294509488e-1, 0.5158423343536001000 ];
        b2 = [ 0.6751786707392436e-1, 0.3215024938520156000, ...
               0.6609491961867980000, 0.2777161669764050000 ];
        w1 = [ 0.2651702815743450e-1, 0.4388140871444811e-1, ...
               0.2877504278497528e-1, 0.6749318700980879e-1 ]*2;  
       
%**************************************************************************
% p = 8                                              (Lower bound = 15 pts)
%**************************************************************************
    elseif p == 8       
    
% 16 Points, PI (Dunavant, verfied)
%--------------------------------------------------------------------------    
        m  = [ 1, 3, 3, 3, 6 ];
        b1 = [ 0.333333333333333, 0.081414823414554, 0.658861384496480, ...
               0.898905543365938, 0.008394777409958 ];
        b2 = [ 0.333333333333333, 0.459292588292723, 0.170569307751760, ...
               0.050547228317031, 0.263112829634638 ];
        w1 = [ 0.144315607677787, 0.095091634267285, 0.103217370534718, ...
               0.032458497623198, 0.027230314174435 ];  
        
% 15 Points, PO (WAY OUT), Optimal (Cools)
%--------------------------------------------------------------------------        
%         m  = [ 3, 3, 3, 3, 3 ];
%         b1 = [ 0.205611832045435,  0.056127355009319, 0.034746808827471,...
%                0.064732904977498, -2.968960232737531 ];
%         b2 = [ 0.281041247315110,  0.630621434318956, 0.313477887523733,...
%                0.870165101563563,  3.623168221569262 ];
%         w1 = [ 0.133881535279832,  0.087819113582442, 0.058571435280332,...
%                0.053061248869561   0.000000000321167 ];    
%**************************************************************************
% p = 9                                              (Lower bound = 17 pts)
%**************************************************************************
    elseif p == 9  
    
% 19 Points, PI (Dunavant)
%--------------------------------------------------------------------------
        m  = [ 1, 3, 3, 3, 3, 6 ];
        b1 = [ 0.3333333333333333333333333333333333,...  
               0.4896825191987376277837069248361928,...
               0.0447295133944527098651065899662764,...
               0.4370895914929366372699303644353550,...
               0.1882035356190327302409612804673356,...        
               0.2219629891607656956751025276931911 ];
        b2 = [ 0.3333333333333333333333333333333333,...  
               0.4896825191987376277837069248361928,...
               0.0447295133944527098651065899662764,...
               0.4370895914929366372699303644353550,...
               0.1882035356190327302409612804673356,...
               0.7411985987844980206900798735234238 ];
        w1 = [ 0.0971357962827988338192419825072886,... 
               0.0313347002271390705368548312872093,...
               0.0255776756586980312616787985589998,...
               0.0778275410047742793167393562994040,...
               0.0796477389272102530328917742640453,...
               0.0432835393772893772893772893772894 ];       
%**************************************************************************
% p = 10                                             (Lower bound = 21 pts)
%**************************************************************************
    elseif p == 10

% 24 Points, PI (Xiao)
%--------------------------------------------------------------------------     
    
% 25 Points, PI (Dunavant)
%--------------------------------------------------------------------------    
    
            m  = [ 1,3,3,3,3,6,6 ];
            b1 = [ 0.3333333333333333333333333333333333,...  
                   0.4272731788467755380904427175154472,...
                   0.1830992224486750205215743848502200,...
                   0.4904340197011305874539712223768484,...
                   0.0125724455515805327313290850210413,...
                   0.6542686679200661406665700955876279,...               
                   0.1228045770685592734301298174812812 ];
            b2 = [ 0.3333333333333333333333333333333333,...  
                   0.4272731788467755380904427175154472,...
                   0.1830992224486750205215743848502200,...
                   0.4904340197011305874539712223768484,...
                   0.0125724455515805327313290850210413,...
                   0.3080460016852477000000000000000000,...
                   0.0333718337393047862408164417747804 ];
            w1 = [ 0.0809374287976228802571131238165019,...
                   0.0772985880029631216825069823803434,...
                   0.0784576386123717313680939208343967,...
                   0.0174691679959294869176071632906781,...
                   0.0042923741848328280304804020901319,...
                   0.0374688582104676429790207654850445,...
                   0.0269493525918799596454494795810967 ];        
    end
%    else
% 22 Points, PO (Cools)
%--------------------------------------------------------------------------
%    end

%**************************************************************************
% p = 11                                             (Lower bound = 24 pts)
%**************************************************************************    
% elseif p == 11    
% % 27 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% % 27 Points, PO (Dunavant)
% %--------------------------------------------------------------------------
% % 28 Points, PI (Lynness)
% %--------------------------------------------------------------------------
% 
% %**************************************************************************
% % p = 12                                             (Lower bound = 28 pts)
% %**************************************************************************
% elseif p == 12
% % 32 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------    
% % 33 Points, PI (Dunavant)
% %--------------------------------------------------------------------------
% %**************************************************************************
% % p = 13                                             (Lower bound = 31 pts)
% %**************************************************************************
% elseif p == 13
% % 35 Points, PI (Xiao, not given in publication)
% %-------------------------------------------------------------------------- 
% % 36 Points, NO (Berntsent)
% %--------------------------------------------------------------------------
% % 37 Points, PI (Dunavant)
% %--------------------------------------------------------------------------
% %**************************************************************************
% % p = 14                                             (Lower bound = 36 pts)
% %**************************************************************************
% elseif p == 14
% % 41 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------     
% % 42 Points, PI (Dunavant)
% %--------------------------------------------------------------------------
% %**************************************************************************
% % p = 15                                             (Lower bound = 40 pts)
% %**************************************************************************
% elseif p == 15
% % 48 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% % 48 Points, PO (Dunavant)
% %--------------------------------------------------------------------------
% %**************************************************************************
% % p = 16                                             (Lower bound = 45 pts)
% %**************************************************************************
% elseif p == 16
% % 52 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% % 52 Points, PO (Dunavant)
% %--------------------------------------------------------------------------
% %**************************************************************************
% % p = 17                                             (Lower bound = 49 pts)
% %**************************************************************************
% elseif p == 17
% % 58 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------    
% % 61 Points, PI (Dunavant)                           
% %--------------------------------------------------------------------------
% %**************************************************************************
% % p = 18                                             (Lower bound = 55 pts)
% %**************************************************************************
% elseif p == 18
% % 65 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------     
% % 70 Points, NO (Dunavant)
% %--------------------------------------------------------------------------
% %**************************************************************************
% % p = 19                                             (Lower bound = 60 pts)
% %**************************************************************************
% elseif p == 19
% % 71 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------     
% % 73 Points, PI (Dunavant)
% %--------------------------------------------------------------------------
% %**************************************************************************
% % p = 20                                             (Lower bound = 66 pts)
% %**************************************************************************
% elseif p == 20
% % 79 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------     
% % 79 Points, NO (Dunavant)
% %--------------------------------------------------------------------------
% % 88 Points, PI (Zhang)
% %--------------------------------------------------------------------------       
% elseif p == 21
% % 85 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------  
% elseif p == 22
% % 95 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 23
% % 102 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 24
% % 110 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 25
% % 118 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 26
% % 127 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 27
% % 136 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 28
% % 147 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 29
% % 156 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 30
% % 168 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 31
% % 181 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 32
% % 193 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 33
% % 204 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 34
% % 214 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 35
% % 228 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 36
% % 243 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 37
% % 252 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 38
% % 267 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 39
% % 282 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 40
% % 295 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 41
% % 309 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 42
% % 324 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 43
% % 339 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 44
% % 354 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 45
% % 370 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 46
% % 385 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 47
% % 399 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 48
% % 423 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 49
% % 435 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
% elseif p == 50
% % 453 Points, PI (Xiao, not given in publication)
% %--------------------------------------------------------------------------
%**************************************************************************
% Exceptions -- Use double product formula
%**************************************************************************    
else
    N = ceil((p+1)/2);
    [X_x,W_x] = gauss_jacobi_quad(N,0,0);
    [X_y,W_y] = gauss_jacobi_quad(N,1,0);
    k = 1;
    for i = 1:N
        for j = 1:N
            X(k) = X_x(i);
            Y(k) = X_y(j);
            w1(k) = W_x(i)*W_y(j)/2;
            m(k) = 1;
            k = k + 1;                
        end
    end
    
    [X',Y']
    % Transform points from the square to the triangle
    X = 1/2*(1+X).*(1-Y) - 1;
    B = cart2bary([X',Y'],[ -1,-1; 1,-1; -1, 1 ]);
    b1 = B(:,1); b2 = B(:,2); b3 = B(:,3); 
end







k = 1;
ng = sum(m);
x = zeros(ng,1);
y = zeros(ng,1);
w = zeros(ng,1);
b3 = 1 - (b1 + b2);
for i = 1:length(m)    
    if m(i) == 1        
        x(k) = a(1)*b1(i) + b(1)*b2(i) + c(1)*b3(i);
        y(k) = a(2)*b1(i) + b(2)*b2(i) + c(2)*b3(i);   
        w(k:k+m(i)-1) = area*w1(i);        
        k = k+1;        
    end   
    if m(i)==2 % Reflection symmetries
        x(k)   = a(1)*b1(i) + b(1)*b2(i) + c(1)*b3(i);
        y(k)   = a(2)*b1(i) + b(2)*b2(i) + c(2)*b3(i);        
        x(k+1) = a(1)*b2(i) + b(1)*b3(i) + c(1)*b1(i);
        y(k+1) = a(2)*b2(i) + b(2)*b3(i) + c(2)*b1(i);               
        w(k:k+m(i)-1) = area*w1(i);        
        k = k+2;                
    end
    if m(i)==3 % Rotation symmetries       
        x(k)   = a(1)*b1(i) + b(1)*b2(i) + c(1)*b3(i);
        y(k)   = a(2)*b1(i) + b(2)*b2(i) + c(2)*b3(i);        
        x(k+1) = a(1)*b2(i) + b(1)*b3(i) + c(1)*b1(i);
        y(k+1) = a(2)*b2(i) + b(2)*b3(i) + c(2)*b1(i);        
        x(k+2) = a(1)*b3(i) + b(1)*b1(i) + c(1)*b2(i);
        y(k+2) = a(2)*b3(i) + b(2)*b1(i) + c(2)*b2(i);        
        w(k:k+m(i)-1) = area*w1(i);        
        k = k+3;        
    end  
    if m(i)==6        
        x(k)   = a(1)*b1(i) + b(1)*b2(i) + c(1)*b3(i);
        y(k)   = a(2)*b1(i) + b(2)*b2(i) + c(2)*b3(i);        
        x(k+1) = a(1)*b2(i) + b(1)*b3(i) + c(1)*b1(i);
        y(k+1) = a(2)*b2(i) + b(2)*b3(i) + c(2)*b1(i);        
        x(k+2) = a(1)*b3(i) + b(1)*b1(i) + c(1)*b2(i);
        y(k+2) = a(2)*b3(i) + b(2)*b1(i) + c(2)*b2(i);        
        x(k+3) = a(1)*b1(i) + b(1)*b3(i) + c(1)*b2(i);
        y(k+3) = a(2)*b1(i) + b(2)*b3(i) + c(2)*b2(i);        
        x(k+4) = a(1)*b2(i) + b(1)*b1(i) + c(1)*b3(i);
        y(k+4) = a(2)*b2(i) + b(2)*b1(i) + c(2)*b3(i);        
        x(k+5) = a(1)*b3(i) + b(1)*b2(i) + c(1)*b1(i);
        y(k+5) = a(2)*b3(i) + b(2)*b2(i) + c(2)*b1(i);        
        w(k:k+m(i)-1) = area*w1(i);        
        k = k+6;
    end
end


% close all
% hold on
% plot([a(1),b(1)],[a(2),b(2)],'k')
% plot([b(1),c(1)],[b(2),c(2)],'k')
% plot([c(1),a(1)],[c(2),a(2)],'k')
% %plot([a(1),-1/3],[a(2),-1/3],'k--')
% %plot([b(1),-1/3],[b(2),-1/3],'k--')
% %plot([c(1),-1/3],[c(2),-1/3],'k--')
% %plot([0,0],[0,1.5],'k:')
% %plot([0,0],[1.5,0],'k')
% 
% plot([a(1),1/2*(b(1)+c(1))],[a(2),1/2*(b(2)+c(2))],'k:')
% plot([b(1),1/2*(a(1)+c(1))],[b(2),1/2*(a(2)+c(2))],'k:')
% plot([c(1),1/2*(a(1)+b(1))],[c(2),1/2*(a(2)+b(2))],'k:')
% 
% for k=1:length(x)
%     plot(x(k),y(k),'bo','MarkerEdgeColor',[ 0.702 0.000 0.000 ],'MarkerFaceColor',[ 0.702 0.000 0.000 ],'MarkerSize',7.0)
% %    text(x(k)+0.05,y(k),int2str(k))
% end
% 
% %axis([-1.0,1.0,-0.25,sqrt(3)])
% axis off
% 
% set(gca,'XTick',[-1.0 0.0 1.0])
% set(gca,'YTick',[-1.0 0.0 1.0])





    
    
    
    
