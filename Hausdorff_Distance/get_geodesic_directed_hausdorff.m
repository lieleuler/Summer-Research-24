% calculate the directed hausdorff distance from geodesic_1 to geodesic_2
function geodesic_directed_hausdorff = get_geodesic_directed_hausdorff(geodesic_1,geodesic_2)
    % mobius transformation for geodesic_1
    [a, b, c, d] = geodesic_1.find_flt_to_imag_axis();
    geodesic_1 = geodesic_1.fractional_linear_transform(a, b, c, d); % transform geodesic_1 to a straight line on imaginary axis
    geodesic_2 = geodesic_2.fractional_linear_transform(a, b, c, d);

    c = geodesic_2.get_center_on_real_line(); % calculate the center of the circle on which geodesic_2 lies
    r = geodesic_2.get_radius_from_center(); % calculate the radius of the circle on which geodesic_2 lies

    syms si % angle at center of d2 from x-axis to the chosen point p on s

    % calculate the range of si based on end points of the segement
    [start_pt_1, end_pt_1] = geodesic_1.get_endpoints();
    [start_pt_2, end_pt_2] = geodesic_2.get_endpoints();
    
    if c >= 0 
        si_1 = pi - atan(imag(start_pt_1)/c);
        si_2 = pi - atan(imag(end_pt_1)/c);
    else
        si_1 = atan(imag(start_pt_1)/-c);
        si_2 = atan(imag(end_pt_1)/-c);
    end
    range_si = sort([si_1, si_2], "ComparisonMethod", "abs");

    % calculate the range of theta based on end points of the segement
    theta_1 = mod(atan((real(start_pt_2)-c)/imag(start_pt_2)), pi);
    theta_2 = mod(atan((real(end_pt_2)-c)/imag(end_pt_2)), pi);
    range_theta = sort([theta_1,theta_2]);

    % the directed hausdorff distance is obtained at the boundaries of geodesic_1
    x_si_1 = (2*r*c*cos(range_si(1))^2) / (c^2*cos(range_si(1))^2 + r^2); % evaluate x = cos(theta) at boundaries of geodesic)_1
    x_si_2 = (2*r*c*cos(range_si(2))^2) / (c^2*cos(range_si(2))^2 + r^2); % evaluate x = cos(theta) at boundaries of geodesic)_2
    theta_si_1 = acos(x_si_1); % evaluate theta at boundaries of geodesic)_1
    theta_si_2 = acos(x_si_2); % evaluate theta at boundaries of geodesic)_2
    
    % calculate the directed hausdorff distance from the start point of geodesic_1
    if theta_si_1 < range_theta(1) || theta_si_1 > range_theta(2) % theta_si_1 attained when outside range_theta
        dist_1 = max(dist_H(start_pt_1,start_pt_2),dist_H(start_pt_1,end_pt_2)); % take boundary points of geodesic_2 instead
    else
        dist_1 = acosh(abs(-1+(2*(c^2 - 2*c*r*x_si_1) * (cos(range_si(1)))^2 + 2*r^2)/(c*r * sqrt(1-x_si_1^2) * sin(2*range_si(1)))));
    end

    % calculate the directed hausdorff distance from the end point boundary of geodesic_1
    if theta_si_2 < range_theta(1) || theta_si_2 > range_theta(2) % theta_si_2 attained when outside range_theta
        dist_2 = max(dist_H(end_pt_1,start_pt_2),dist_H(end_pt_1,end_pt_2)); % take boundary points of geodesic_2 instead
    else
        dist_2 = acosh(abs(-1+(2*(c^2 - 2*c*r*x_si_2) * (cos(range_si(2)))^2 + 2*r^2)/(c*r * sqrt(1-x_si_1^2) * sin(2*range_si(2)))));
    end
    
    % output the directed hausdorff distance
    geodesic_directed_hausdorff = max(dist_1,dist_2);
end
