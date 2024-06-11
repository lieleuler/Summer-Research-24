
function circleEq = get_circle_eq(center, radius)   
    % Construct symbolic equation of the circle
    circleEq = @(s, t) dist_H(s + t*i,center) == radius;
end
