% calculate the angle between intersecting points
function rangeTn = getRangeTn(t, intersection, segment, p_radius, seg_radius)
    % Filter down to 1-2 intersections
    if height(intersection) > 2
        "HEYYYYY"
        intersection
        new_intersection = zeros(0, 1);
        for int = intersection.'
            assisstance_line = GeodesicSegment(t, int);
            if ~assisstance_line.intersects_geodesic(segment, false, true)
                new_intersection = [new_intersection; int];
            end
        end
        intersection = new_intersection;
    end
    if height(intersection) ~= 2
        [e1, e2] = segment.get_endpoints();
        [t, e1, e2, p_radius, seg_radius]
        error("More than two valid intersections detected, even " + ...
            "after filitering (" + height(intersection) + " ints)")
    end

    % No
    u1 = intersection(1); 
    u2 = intersection(2);
    
    % calculate the slope of the tangent of a circle
    slope_1 = getTanSlope(u1,t);
    slope_2 = getTanSlope(u2,t);
    
    % calculate the angle by setting vertical line as angle zero
    t_x = real(t);
    u1_x = real(u1);
    u2_x = real(u2);

    % TO-DO: Case with slope == 0
    if slope_1 > 0 && u1_x > t_x
        theta_1 = 3*pi/2 + atan(slope_1); % first quadrant
    elseif slope_1 < 0 && u1_x <= t_x
        theta_1 = pi/2 + atan(slope_1); % second quadrant  
    elseif slope_1 > 0 && u1_x <= t_x
        theta_1 = pi/2 + atan(slope_1); % third quadrant
    else
        theta_1 = 3*pi/2 + atan(slope_1); % fourth quadrant
    end

    if slope_2 > 0 && u2_x > t_x
        theta_2 = 3*pi/2 + atan(slope_2); % first quadrant
    elseif slope_2 < 0 && u2_x <= t_x
        theta_2 = pi/2 + atan(slope_2); % second quadrant  
    elseif slope_2 > 0 && u2_x <= t_x
        theta_2 = pi/2 + atan(slope_2); % third quadrant
    else
        theta_2 = 3*pi/2 + atan(slope_2); % fourth quadrant
    end
    
    %if theta_1 == theta_2
        %"HHHH"
        %[e1, e2] = segment.get_endpoints();
        %data = [t, e1, e2, p_radius, seg_radius]
        %intersection
    %end
    
    % Calculate the angle range outside of boundary
    if theta_2 <= theta_1
        theta_2_adjusted = theta_2 + 2*pi;
    else
        theta_2_adjusted = theta_2;
    end
    test_angle = (theta_1 + theta_2_adjusted)/2;
    test_point = get_point_along_direction(t, test_angle, p_radius);
    % CWW from theta_1
    if segment.dist_from_point(test_point) >= seg_radius
        rangeTn = [theta_1, theta_2_adjusted];
    % CWW from theta_2
    else
        if theta_1 <= theta_2
            theta_1_adjusted = theta_1 + 2*pi;
        else
            theta_1_adjusted = theta_1;
        end
        rangeTn = [theta_2, theta_1_adjusted];
    end
end

% calculate the slope of the tangent of a circle
function tanSlope = getTanSlope(point,center)
    u = point;
    t = center;
    c = getGeoCenter(u,t); % center of the great geodesic circle

    % Evaluate the derivative of circle with center c and point t
    % at the center of step size circle t
    tanSlope = (real(t) - real(c)) / (imag(c) - imag(t));
end
