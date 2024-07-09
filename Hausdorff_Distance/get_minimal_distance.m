

% calculate the minimal distance MD(geo_1,geo_2) for pruning
function minimal_distance = get_minimal_distance(geodesic_1,geodesic_2)    
    % mobius transformation for geodesic_1
    [a, b, c, d] = geodesic_1.find_flt_to_imag_axis();
    geodesic_1 = geodesic_1.fractional_linear_transform(a, b, c, d); % transform geodesic_1 to a straight line on imaginary axis
    geodesic_2 = geodesic_2.fractional_linear_transform(a, b, c, d);
    
    c = geodesic_2.get_center_on_real_line() % calculate the center of the circle on which geodesic_2 lies
    r = geodesic_2.get_radius_from_center() % calculate the radius of the circle on which geodesic_2 lies
        
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

    % Determine si value in [0, pi] which minimizes the distance between
    % the two geodesics, and save it as si_root
    si_root = -atan(sqrt(abs(r^2 - c^2))/c);
    if c > 0
        si_root = pi - si_root;
    end

    % case 1: si_root is outside angles bounds range_si
    if si_root < range_si(1) || range_si(2) < si_root
        % calculate projection of nodes of geodesic_2
        [node, dismissal] = sort([start_pt_2,end_pt_2],"ComparisonMethod","real");
        node_projection = sqrt(imag(node)^2 + real(node)^2) * 1i;
 
        [lower_pt_1,higher_pt_1] = sort([start_pt_1, end_pt_1],"ComparisonMethod","abs");

        % case 1a
        if imag(node_projection) < imag(lower_pt_1) || imag(higher_pt_1) < imag(node_projection)
            minimal_distance = min(dist_H(start_pt_1,node), dist_H(end_pt_1,node));

        % case 1b
        else
            slope_orth_node = -1/get_tan_slope(node,c);
            p_node = imag(node) - slope_orth_node * real(node); % where start_pt_2 is the projection of node
            minimal_distance = dist_H(p_node,node);
        end
    
    % Calculate the minimizing
    x_si_root = (-2*c*r*cos(si_root)^2) / (c^2 + r^2*cos(si_root)^2);
    theta_si = acos(x_si_root); % evaluate numerical value of theta at si_root

    % case 2: si_root is within angle bound range_si BUT theta_si is obtained outside of the angle range range_theta 
    elseif  theta_si > range_theta(2) || theta_si < range_theta(1)
        % calculate projection of nodes of geodesic_2
        [node, dismissal] = sort([start_pt_2,end_pt_2],"ComparisonMethod","real");
        node_projection = sqrt(imag(node)^2 + real(node)^2) * 1i;
 
        [lower_pt_1,higher_pt_1] = sort([start_pt_1, end_pt_1],"ComparisonMethod","abs");
        
        % case 2b
        if imag(node_projection) < imag(lower_pt_1) || imag(higher_pt_1) < imag(node_projection)
            minimal_distance = min(dist_H(start_pt_1,node), dist_H(end_pt_1,node));

        % case 2a
        else
            minimal_distance = dist_H(node_projection,node);
        end

    % case 3: both in
    else 
        numerator_at_root = 2*c^2 + 2*r^2*cos(si_root)^2 + 4*c*r*cos(si_root)*(x_si_root * cos(si_root) + sqrt(1 - x_si_root^2)*sin(si_root));
        denominator_at_root = 2 * c*r*sin(2*si_root) * sqrt(1- x_si_root^2);
        minimal_distance = acosh(1 - numerator_at_root/denominator_at_root);
    end
end

function min_dist = min_distance_from_si(si, c, r)
    x_si = (-2*c*r*cos(si)^2) / (c^2 + r^2*cos(si)^2); % evaluate numerical value of x = cos(theta) at si_root
    numerator_at_root = 2*c^2 + 2*r^2*cos(si)^2 + 4*c*r*cos(si_root)*(x_si * cos(si) + sqrt(1 - x_si^2)*sin(si));
    denominator_at_root = 2 * c*r*sin(2*si) * sqrt(1- x_si^2);
    min_dist = acosh(1 - numerator_at_root/denominator_at_root);
end