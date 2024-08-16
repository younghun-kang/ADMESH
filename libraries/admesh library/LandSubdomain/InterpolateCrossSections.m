function CS = InterpolateCrossSections(xyzFun,LineXY)

for i = 1 : size(LineXY,1)       
    px = LineXY(i,1);
    py = LineXY(i,2);
    if i == 1
        tx = LineXY(i+1,1) - LineXY(i,1);
        ty = LineXY(i+1,2) - LineXY(i,2);
    elseif i == size(LineXY,1)
        tx = LineXY(i,1) - LineXY(i-1,1);
        ty = LineXY(i,2) - LineXY(i-1,2);
    else
        tx = LineXY(i+1,1) - LineXY(i-1,1);
        ty = LineXY(i+1,2) - LineXY(i-1,2);
    end
    nx = ty/sqrt(tx^2 + ty^2);
    ny = -tx/sqrt(tx^2 + ty^2);
    
    % Linear interpolation.    
    fx = @(d) px + nx*d;
    fy = @(d) py + ny*d;
    fz = @(d) xyzFun(px+nx*d,py+ny*d);

    CS(i).fx = fx;
    CS(i).fy = fy;
    CS(i).fz = fz;    
end

end







