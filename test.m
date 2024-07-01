
% 0         6.2832
% 5.5523    5.7252
% 5.3274    5.4668
% 5.1664    5.5308
% 10.7724   12.4912


e1 = 0.0000 + 1.0000i;
e2 = 0.0242 + 1.1022i;
%e1 = 4.8859i;
%e2 = 4.4210i;
g = GeodesicSegment(e1, e2);

z = 0.1130 + 1.1733i;

ranges = zeros(0, 2);
for i = 0:10
    t_i = 0.1/11 * i;
    t_i_next = 0.1/11 * (i + 1);
    sub_g = GeodesicSegment(g.travel_from_start(t_i), g.travel_from_start(t_i_next));

    lowerBd = (1/1.1) * (4 - (1 + i/11)) * 0.1;
    fat = (1 - 1/1.1) * 0.1;
    s = lowerBd + fat;

    intersection_i = intersections_of_point_and_segment_ngbhs(z, sub_g, ...
                0.1, s);
    inside = (sub_g.dist_from_point(z) <= s);
    range_i = getRangeTn(z, intersection_i, inside, sub_g);
    ranges = [ranges;range_i];
end

%merged_range = MergeRange(ranges);


[a, b, c, d] = g.find_flt_to_imag_axis();
[x, w] = g.fractional_linear_transform(a, b, c, d).get_endpoints()
y = (a*z + b) / (c*z + d)
dist_H(w, y);


%h = dist_H(0.0000 + 1.0000i,  0.1396 + 0.8002i);
%d = 3/1.1 * 0.1 + (1 - 1/1.1) * 0.1;
%h = dist_H(z, g.travel_from_start(0.08));
%lower_bound = (4 - 1)/1.1 * 0.1;

%e_sub_1 = g.travel_from_start(0.1/11 * 1);
%e_sub_2 = g.travel_from_start(0.1/11 * 2);

%g_sub = GeodesicSegment(e_sub_1, e_sub_2);

%int = intersections_of_point_and_segment_ngbhs(z, g_sub, 0.1, 0.1 * (4 - (1 + 0.1/11))/1.1 + (1 - 1/1.1)*0.1);

%dist_H(e2, z);
%getRangeTn(z, int);


%nudge = 0;
%l = g.get_length
%dist_H(z, g.travel_from_start(nudge))
%bound = (0.1*(3 - 1) + nudge)/1.1

%dist_H(e1, -0.0114 +0.90549i)

%get_point_along_direction(z, 1.9412 + pi/2, 0.1)

%getTanSlope(i2, z)
%getRangeTn(z, [i1, i2])

%dist_H(z, i2)
%g.dist_from_point(i1)

%double(intersections_of_point_and_segment_ngbhs(z, GeodesicSegment(e1, e2), 0.1, 0.2727))

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

    u = point;
    t = center;
    c = getGeoCenter(u,t); % center of the great geodesic circle

    % Evaluate the derivative of circle with center c and point t
    % at the center of step size circle t
    tanSlope = (real(t) - real(c)) / (imag(c) - imag(t));
end


function rangeTn = getRangeTn(t, intersection, t_inside_boundary, segment)
    % Filter down to 1-2 intersections
    if isempty(intersection)
        rangeTn = [0, 2*pi];
        return
    elseif height(intersection) > 2
        "HEYYYYY"
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
            "after filitering")
    end

    u1 = intersection(1); 
    u2 = intersection(2); % TO-DO: Case for 4 intersections
    
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


% merge multiple ranges
function mergedRange = MergeRange(ranges)
    % Initialize the intersection range
    mergedRange = [-inf, inf];
    
    % Find the intersection of all ranges
    for i = 1:size(ranges, 1)
        curr_range = ranges(i, :);
        if curr_range(1) > 2*pi
            curr_range = curr_range - [2*pi, 2*pi];
        end
        if curr_range(2) > 2*pi
            curr_range = [curr_range(1), 2*pi;
                          0, curr_range(2) - 2*pi];
        end
        mergedRange_size = size(mergedRange, 1);
        curr_range_size = size(curr_range, 1);
        new_mergedRange = zeros(0, 2);
        for j = 1:mergedRange_size
            for k = 1:curr_range_size
                r_lower = max(mergedRange(j, 1), curr_range(k, 1));
                r_upper = min(mergedRange(j, 2), curr_range(k, 2));
                if r_upper > r_lower
                    new_mergedRange = [new_mergedRange ; [r_lower, r_upper]];
                end
            end
        end
        mergedRange = new_mergedRange;
    end
    
    % Check if the intersection is valid
    if isempty(mergedRange)
        double(ranges)
        error('The ranges do not overlap.');
    end
end

% generate bounded t_n using randomization
function tn = generateTn(t,range,step_size) % ?might be able to generate "geodesic"
    
    % Generate random value within size of range for uniform distribution
    range_size = 0;
    for i = 1:size(range, 1)
        range_size = range_size + (range(i, 2) - range(i, 1));
    end
    random_value = range_size * rand();

    % Get phi via random value
    for i = 1:size(range, 1)
        random_value = random_value - (range(i, 2) - range(i, 1));
        if random_value <= 0
            phi = range(i, 2) + random_value;
            break
        end
    end

    % Calculate the coordinates of the random point
    % NOTE: angles in range place 0 vertically up, while angles in the call
    % to this function place 0 to the right like normal polar coordinates.
    % Thus, we add pi/2 for conversion.
    tn = get_point_along_direction(t, double(phi + pi/2), step_size);
end
