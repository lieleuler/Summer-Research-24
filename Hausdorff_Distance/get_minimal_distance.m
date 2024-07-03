% calculate the minimal distance MD(geo_1,geo_2) for pruning
function minimal_distance = get_minimal_distance(geodesic_1,geodesic_2) 
    % mobius transformation for geodesic_1
    [a, b, c, d] = geodesic_1.find_flt_to_imag_axis();
    geodesic_1 = geodesic_1.fractional_linear_transform(a, b, c, d); % transform geodesic_1 to a straight line on imaginary axis
    geodesic_2 = geodesic_2.fractional_linear_transform(a, b, c, d);
    
    c = geodesic_2.get_center_on_real_line(); % calculate the center of the circle on which geodesic_2 lies
    r = geodesic_2.get_radius_from_center(); % calculate the radius of the circle on which geodesic_2 lies

    syms si % angle at center of d2 from x-axis to the chosen point p on s

    % calculate the range of si based on end points of the segement
    [start_pt_1, end_pt_1] = geodesic_1.get_endpoints();
    
    if c >= 0 
        si_1 = pi - atan(imag(start_pt_1)/c);
        si_2 = pi - atan(imag(end_pt_1)/c);
    else
        si_1 = atan(imag(start_pt_1)/-c);
        si_2 = atan(imag(end_pt_1)/-c);
    end
    range_si = [si_1, si_2];

    % calculate the range of theta based on end points of the segement
    [start_pt_2, end_pt_2] = geodesic_2.get_endpoints();

    theta_1 = mod(atan(imag(start_pt_2)/(real(start_pt_2)-c)), pi);
    theta_2 = mod(atan(imag(end_pt_2)/(real(end_pt_2)-c)), pi);
    range_theta = sort([theta_1,theta_2]);

    x_si = (-2*c*r*cos(si)^2) / (c^2 + r^2*cos(si)^2);
    numerator = 2*c^2 + 2*r^2*cos(si)^2 + 4*c*r*cos(si)*(x_si * cos(si) + sqrt(1 - x_si^2)*sin(si));
    denominator = 2 * c*r*sin(2*si) * sqrt(1- x_si^2);
    g = 1 - numerator/denominator;
    g_prime = diff(g, si);

    si_roots = vpasolve(g_prime == 0, si, sort(range_si));
    
    % output MD
    if isempty(si_roots) % no critical point in the specified range_si
        minimal_distance = min([dist_H(start_pt_1,start_pt_2),dist_H(start_pt_1,end_pt_2),dist_H(start_pt_2,end_pt_1),dist_H(end_pt_1,end_pt_2)]);
        return
    end

    x_si = (-2*c*r*cos(si_roots)^2) / (c^2 + r^2*cos(si_roots)^2); % evaluate numerical value of x = cos(theta) at si_roots
    theta_si = acos(x_si); % evaluate numerical value of theta at si_roots
    % If the minimum is obtained outside of the angle range of geodesic_2
    if theta_si > range_theta(2) || theta_si < range_theta(1) % theta_si outside the range_theta 
        minimal_distance = min([dist_H(start_pt_1,start_pt_2),dist_H(start_pt_1,end_pt_2),dist_H(start_pt_2,end_pt_1),dist_H(end_pt_1,end_pt_2)]);
    % If the minimum is obtained inside of the angle range of geodesic_2
    else % normal case 
        numerator = 2*c^2 + 2*r^2*cos(si_roots)^2 + 4*c*r*cos(si_roots)*(x_si * cos(si_roots) + sqrt(1 - x_si^2)*sin(si_roots));
        denominator = 2 * c*r*sin(2*si_roots) * sqrt(1- x_si^2);
        minimal_distance = acosh(1 - numerator/denominator);
    end
end
