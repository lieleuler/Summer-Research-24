% calculate the slope of the tangent of a circle
function tanSlope = get_tan_slope(point,center)
    syms x

    u = point;
    t = center;
    c = getGeoCenter(u,t); % center of the great geodesic circle
    R = dist_H(u,c); % radius of the great geodesic circle

    % Express the great geodesic circle as a function
    circFunc = getCircleFunc(c, R) 
    
    % Compute the derivative
    df = diff(circFunc, x);
    
    % Evaluate the derivative at the center of step size circle t
    tanSlope = subs(df, x, t);
    
end
