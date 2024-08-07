% calculate the angle between intersecting points
function rangeTn = getRangeTn(t, intersection, segment, p_radius, seg_radius, is_lower_bound)
    % Filter down to 1-2 intersections
    if height(intersection) > 2
        new_intersection = zeros(0, 1);
        for int = intersection.'
            assisstance_line = GeodesicSegment(t, int);
            if ~assisstance_line.intersects_geodesic(segment, false, false)
                new_intersection = [new_intersection; int];
            end
        end
        intersection = new_intersection;
    end
    if height(intersection) ~= 2
        [e1, e2] = segment.get_endpoints();
        disp("t: " + t)
        disp("e1: " + e1)
        disp("e2: " + e2)
        disp("p_radius: " + p_radius)
        disp("seg_radius: " + seg_radius)
        disp(real(intersection))
        disp(imag(intersection))
        disp(segment.dist_from_point(intersection(1)))
        disp(segment.dist_from_point(intersection(3)))
        error("More than two valid intersections detected, even " + ...
            "after filtering (" + height(intersection) + " ints)")
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

    % Calculate the angle range outside of boundary
    if theta_1 <= theta_2
            theta_1_adjusted = theta_1 + 2*pi;
    else
        theta_1_adjusted = theta_1;
    end
    test_angle_CW = (theta_2 + theta_1_adjusted)/2;
    test_point_CW = get_point_along_direction(t, test_angle_CW, p_radius);
    if theta_2 <= theta_1
        theta_2_adjusted = theta_2 + 2*pi;
    else
        theta_2_adjusted = theta_2;
    end
    test_angle_CWW = (theta_1 + theta_2_adjusted)/2;
    test_point_CWW = get_point_along_direction(t, test_angle_CWW, p_radius);
    % CWW from theta_1
    if (is_lower_bound && segment.dist_from_point(test_point_CWW) >= segment.dist_from_point(test_point_CW)) ...
       || (~is_lower_bound && segment.dist_from_point(test_point_CWW) <= segment.dist_from_point(test_point_CW))
        rangeTn = [theta_1, theta_2_adjusted];
    % CWW from theta_2
    else
        rangeTn = [theta_2, theta_1_adjusted];
    end

    % TESTING
    % if abs(rangeTn(1) - rangeTn(2)) <= 1e-5
    %     "HJK"
    %     [e1, e2] = segment.get_endpoints();
    %     disp("t: " + t)
    %     disp("e1: " + e1)
    %     disp("e2: " + e2)
    %     disp("p_radius: " + p_radius)
    %     disp("seg_radius: " + seg_radius)
    % end

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
