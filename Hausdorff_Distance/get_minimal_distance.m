% calculate the minimal distance MD(geo_1,geo_2) for pruning
function minimal_distance = get_minimal_distance(geodesic_1,geodesic_2,c,r) 
    % mobius transformation for geodesic_1
    geodesic_1 = mobius_transform(geodesic_1); % transform geodesic_1 to a straight line on imaginary axis
    
    c = center(geodesic_2); % calculate the center of the circle on which geodesic_2 lies
    r = radius(geodesic_2); % calculate the radius of the circle on which geodesic_2 lies

    syms si % angle at center of d2 from x-axis to the chosen point p on s

    % calculate the range of si based on end points of the segement
    start_pt_1 = get_start_point(geodesic_1);
    end_pt_1 = get_end_point(geodesic_1);
    
    range_si = sort(atan(-start_pt_(2)/c),atan(-end_pt_1(2)/c));

    % calculate the range of theta based on end points of the segement
    start_pt_2 = get_start_point(geodesic_2);
    end_pt_2 = get_end_point(geodesic_2);

    range_theta = sort(atan((start_pt_2(1)-c)/start_pt_2(2)),atan((end_pt_2(1)-c)/end_pt_2(2)));

    % take derivative of distance function
    % Expression for x = cos(theta)
    x = (2*r*c*cos(si)^2) / (c^2*cos(si)^2 + r^2);
    % Derivative of x with respect to si
    dx_dsi = diff(x, si);

    % Substitute x and x' in the given expression
    combined_expr = (1/(2*c*r*(1 - x^2)^(3/2))) * csc(2*si)^2 * (4*c*(c - 2*r*x)*(1 - x^2) ...
        + 4*(c^2 + 2*r^2 - 2*c*r*x)*(1 - x^2)*cos(2*si) ...
        + 2*(-2*c*r + c^2*x + 2*r^2*x)*dx_dsi*sin(2*si) ...
        + c*(-2*r + c*x)*dx_dsi*sin(4*si));

    % Simplify the expression
    combined_expr = simplify(combined_expr);

    % Define the equation to solve
    eqn = combined_expr == 0;

    % Use numerical solver to find si in the range range_si
    si_roots = vpasolve(eqn, si, range_si);
    x_si = (2*r*c*cos(si_roots)^2) / (c^2*cos(si_roots)^2 + r^2); % evaluate numerical value of x = cos(theta) at si_roots
    theta_si = acos(x_si); % evaluate numerical value of theta at si_roots
    
    % output MD
    if isempty(si_roots) % no critical point in the specified range_si
        minimal_distance = min(dist_H(start_pt_1,start_pt_2),dist_H(start_pt_1,end_pt_2),dist_H(start_pt_2,end_pt_1),dist_H(end_pt_1,end_pt_2));
    elseif theta_si > range_theta(2) || theta_si < range_theta(1) % theta_si outside the range_theta 
        minimal_distance = min(dist_H(start_pt_1,start_pt_2),dist_H(start_pt_1,end_pt_2),dist_H(start_pt_2,end_pt_1),dist_H(end_pt_1,end_pt_2));
    else % normal case 
        minimal_distance = acosh(-1+(2*(c^2 - 2*c*r*x_si) * (cos(si_roots))^2 + 2*r^2)/(c*r * sqrt(1-x_si^2) * sin(2*si_roots)));
    end
end
