% calculate the center of the great geodesic circle given two points on it
function geoCenter = getGeoCenter(s,t)
    x1 = real(s);
    y1 = imag(s);
    x2 = real(t);
    y2 = imag(t);
        
    geoCenter = (x1^2 + y1^2 - x2^2 - y2^2) / (2 * (x1 - x2));
end

% calculate the slope of the tangent of a circle
function tanSlope = getTanSlope(point,center)
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


% calculate the angle between intersecting points
function rangeTn = getRangeTn(t,intersection)
    u1 = intersection(1) 
    u2 = intersection(2) % assuming two intersections. possible modification needed!
    
    % calculate the slope of the tangent of a circle
    slope_1 = getTanSlope(u1,t);
    slope_2 = getTanSlope(u2,t);
    
    % calculate the angle by setting vertical line as angle zero
    if u1.x < t.x
        theta_1 = pi/2 + atan(slope_1); % negative slope gives negative atan  
    else
        theta_1 = 3*pi/2 + atan(slope_1); 
    end

    if u2.x < t.x
        theta_2 = pi/2 + atan(slope_2); % negative slope gives negative atan  
    else
        theta_2 = 3*pi/2 + atan(slope_2); 
    end
    
    % calculate the intersecting angle
    angDiff = abs(theta_1 - theta_2);
    if angDiff < pi % then start from larger-angled arm
        rangeTn_start = max(theta_1,theta_2);
        rangeTn = [rangeTn_start, rangeTn_start + angDiff];
    else
        rangeTn_start = min(theta_1,theta_2);
        rangeTn = [rangeTn_start, rangeTn_start + angDiff];


