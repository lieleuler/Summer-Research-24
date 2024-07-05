% calculate the angle between intersecting points
function rangeTn = getRangeTn(t, intersection, t_inside_boundary, segment)
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
    if u1_x < t_x
        theta_1 = pi/2 + atan(slope_1); % negative slope gives negative atan  
    else
        theta_1 = 3*pi/2 + atan(slope_1);
    end

    if u2_x < t_x
        theta_2 = pi/2 + atan(slope_2); % negative slope gives negative atan  
    else
        theta_2 = 3*pi/2 + atan(slope_2); 
    end            
    
    % Calculate the angle range outside of boundary
    angDiff = abs(theta_1 - theta_2);

    int_geod = GeodesicSegment(intersection(1), intersection(2));
    segment_endpoints = segment.get_endpoints();
    useSmallerAngleDiff = t_inside_boundary ... 
        && sign(imag(intersection(1))) == sign(imag(intersection(2)))  ...
        && ~int_geod.intersects_geodesic(GeodesicSegment(t, segment_endpoints(1)), false, false);

    if useSmallerAngleDiff
        if angDiff < pi % then start from larger-angled arm
            rangeTn_start = min(theta_1,theta_2);
            rangeTn = [rangeTn_start, rangeTn_start + angDiff];
        else
            rangeTn_start = max(theta_1,theta_2);
            rangeTn = [rangeTn_start, rangeTn_start + 2*pi - angDiff];
        end
    else 
        if angDiff < pi % then start from larger-angled arm
            rangeTn_start = max(theta_1,theta_2);
            rangeTn = [rangeTn_start, rangeTn_start + 2*pi - angDiff];
        else
            rangeTn_start = min(theta_1,theta_2);
            rangeTn = [rangeTn_start, rangeTn_start + angDiff];
        end
    end
end
