function IN = PointsInDomain2(X,Y,PTS)

% Loop over each polygon

% Test first polygon
IN = InPolygon(X,Y,PTS.Poly(1).x,PTS.Poly(1).y);

% Test Remaining polygons
for k = 2:length(PTS.Poly)
    
    in = InPolygon(X,Y,PTS.Poly(k).x,PTS.Poly(k).y);
    
    % If (x,y) are in an interior polygon we must remove it
    IN(in) = 0;
    
end

% Younghun: bankline constraints
for k = 1:length(PTS.Constraints)
    if PTS.Constraints(k).num == 18
        if isequal(PTS.Constraints(k).xy(1,:),PTS.Constraints(k).xy(end,:))
            in = InPolygon(X,Y,PTS.Constraints(k).xy(:,1),PTS.Constraints(k).xy(:,2));
            
            % If (x,y) are in an interior polygon we must remove it
            IN(in) = 0;
        end
    end
    
end

end