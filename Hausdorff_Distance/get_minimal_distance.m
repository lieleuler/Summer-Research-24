

% calculate the minimal distance MD(geo_1,geo_2) for pruning
function minimal_distance = get_minimal_distance(geodesic_1,geodesic_2)    
    % mobius transformation for geodesic_1
    [a, b, C, d] = geodesic_1.find_flt_to_imag_axis();
    geodesic_1 = geodesic_1.fractional_linear_transform(a, b, C, d); % transform geodesic_1 to a straight line on imaginary axis
    geodesic_2 = geodesic_2.fractional_linear_transform(a, b, C, d);
    
    c = geodesic_2.get_center_on_real_line(); % calculate the center of the circle on which geodesic_2 lies
    r = geodesic_2.get_radius_from_center(); % calculate the radius of the circle on which geodesic_2 lies
        
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
    range_si_sorted = sort(range_si);
    % disp("range_si: "+range_si_sorted);
    
   % calculate the range of theta based on end points of the segement
    [start_pt_2, end_pt_2] = geodesic_2.get_endpoints();
    diff_start = real(start_pt_2)-c;
    diff_end = real(end_pt_2)-c;
    
    if diff_start >= 0 
        theta_1 = atan(imag(start_pt_2)/diff_start);
        
    else
        theta_1 = pi - atan(imag(start_pt_2)/abs(diff_start));
    end
    
    if diff_end >= 0
        theta_2 = atan(imag(end_pt_2)/diff_end);
    else
        theta_2 = pi - atan(imag(start_pt_2)/abs(diff_end));
    end
    
    range_theta = [theta_1,theta_2];
    range_theta_sorted = sort(range_theta);
    % disp("range_theta: "+range_theta_sorted);

    % determine si value in [0, pi] which minimizes the distance between
    % the two geodesics, and save it as si_root
    si_root = -atan(sqrt(abs(r^2 - c^2))/c);
    if c > 0
        si_root = pi - si_root;
    end

    % calculate the minimizing
    x_si_root = (-2*c*r*cos(si_root)^2) / (c^2 + r^2*cos(si_root)^2);
    theta_si = acos(x_si_root); % evaluate numerical value of theta at si_root
    
    % case 0: intersect
    if geodesic_1.intersects_geodesic(geodesic_2,false,false)
        % disp("under case 0:");
        minimal_distance = 0;
        return


    % case 1: si_root is outside angles bounds range_si
    elseif (si_root < range_si_sorted(1)) || (range_si_sorted(2) < si_root)
        % disp("under case 1: si_root outside range_si");
        
        % choose the end point of interest of geodesic_1 (i.e. tip)
        if abs(si_root - range_si(1)) < abs(si_root - range_si(2)) % unsorted 1=start, 2=end
            tip = start_pt_1;
        else
            tip = end_pt_1;
        end

        % projection of tip to arc
        [a,b,C,d] = geodesic_2.find_flt_to_imag_axis();
        tip_transformed = (a*tip + b)/(C*tip + d);
        tip_projection_transformed = sqrt(real(tip_transformed)^2 + imag(tip_transformed)^2) * 1i;
        tip_projection = (d*tip_projection_transformed - b)/(-C*tip_projection_transformed + a);
        % disp("tip_projection: "+tip_projection);
        
        diff = c - real(tip_projection);
        if diff > 0
            theta_tip_projection = pi - atan(imag(tip_projection)/abs(diff));
        else
            theta_tip_projection = atan(imag(tip_projection)/abs(diff));
        end
        
        % choose the end point of interest of geodesic_2 (i.e. node)
        if abs(theta_tip_projection - range_theta(1)) < abs(theta_tip_projection - range_theta(2)) % unsorted 1=start, 2=end
            node = start_pt_2;
        else
            node = end_pt_2;
        end
        
        % projection of nodes to axis
        node_projection = sqrt(real(node)^2 + imag(node)^2) * 1i;
        % disp("node_projection: "+node_projection);
        if c > 0
            si_node_projection = pi - atan(imag(node_projection)/c);
        else
            si_node_projection = atan(imag(node_projection)/abs(c));
        end

        
        % check if tip_projection falls in geodesic segment
        if (range_theta_sorted(1) < theta_tip_projection) && (theta_tip_projection < range_theta_sorted(2))
            % disp("tip in");
            % disp("theta_tip_projection: "+theta_tip_projection);
            minimal_distance = dist_H(tip,tip_projection);
        else
            % disp("tip out");
            % check if node_projection falls in geodesic segment
            if (range_si_sorted(1) < si_node_projection) && (si_node_projection < range_theta_sorted(2))
                % disp("node in");
                minimal_distance = dist_H(node,node_projection);
            else
                % disp("node out");
                minimal_distance = dist_H(tip,node);
            end
        end
            
    
        
    % case 2: si_root is within angle bound range_si BUT theta_si is obtained outside of the angle range range_theta 
    elseif  (theta_si > range_theta_sorted(2)) || (theta_si < range_theta_sorted(1))
        % disp("under case 2: si_root inside range_si but theta_si outside range_theta");
        % choose the end point of interest of geodesic_1 (i.e. tip)
        if abs(si_root - range_si(1)) < abs(si_root - range_si(2)) % unsorted 1=start, 2=end
            tip = start_pt_1;
        else
            tip = end_pt_1;
        end

        % projection of tip to arc
        [a,b,C,d] = geodesic_2.find_flt_to_imag_axis();
        tip_transformed = (a*tip + b)/(C*tip + d);
        tip_projection_transformed = sqrt(real(tip_transformed)^2 + imag(tip_transformed)^2) * 1i;
        tip_projection = (d*tip_projection_transformed - b)/(-C*tip_projection_transformed + a);

        diff = c - real(tip_projection);
        if diff > 0
            theta_tip_projection = pi - atan(imag(tip_projection)/abs(diff));
        else
            theta_tip_projection = atan(imag(tip_projection)/abs(diff));
        end

        % choose the end point of interest of geodesic_2 (i.e. node)
        if abs(theta_tip_projection - range_theta(1)) < abs(theta_tip_projection - range_theta(2)) % unsorted 1=start, 2=end
            node = start_pt_2;
        else
            node = end_pt_2;
        end


        % projection of nodes of geodesic_2 to axis
        if abs(theta_tip_projection - range_theta(1)) < abs(theta_tip_projection - range_theta(2)) % unsorted 1=start, 2=end
            node = start_pt_2;
        else
            node = end_pt_2;
        end
        
        node_projection = sqrt(real(node)^2 + imag(node)^2) * 1i;
        if c > 0
            si_node_projection = pi - atan(imag(node_projection)/c);
        else
            si_node_projection = atan(imag(node_projection)/abs(c));
        end

        % since theta_si outside range_theta, tip_projection must be outside range_theta, so no discussion needed
        % only check if tip_projection falls in geodesic segment 
        if (range_si_sorted(1) < si_node_projection) && (si_node_projection < range_theta_sorted(2))
            % disp("node in");
            minimal_distance = dist_H(node,node_projection);
        else
            % disp("node out");
            minimal_distance = dist_H(tip,node);
        end


    % case 3: both in
    else 
        % disp("under case 3: si_root inside range_si and theta_si inside range_theta");
        % disp("si root: "+si_root);
        numerator_at_root = 2*c^2 + 2*r^2*cos(si_root)^2 + 4*c*r*cos(si_root)*(x_si_root * cos(si_root) + sqrt(1 - x_si_root^2)*sin(si_root));
        denominator_at_root = 2 * c*r*sin(2*si_root) * sqrt(1- x_si_root^2);
        minimal_distance = acosh(1 - numerator_at_root/denominator_at_root);
    end
end

