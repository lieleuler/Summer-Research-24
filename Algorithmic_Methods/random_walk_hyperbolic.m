
function [points, ranges, phi] = random_walk_hyperbolic(n_steps, lambda, eps, ...
    step_size, min_segment_splits)
    % This function performs a random walk in the hyperbolic upper half-plane
    % n_steps: Number of steps in the random walk
    % step_size: Fixed hyperbolic distance to move at each step

    % Set up parameters
    validate_parameter_conditions(lambda, eps, step_size, n_steps, min_segment_splits)

    % Initialize the starting point
    z = 1i;  % Start at (0, 1) to avoid the real axis
    
    % Pre-allocate array for plotting
    points = zeros(n_steps + 1, 1);
    points(1) = z;

    % Pre-calculate the smallest angle allowed by the lower bound of the
    % segment adjacent to the point 
    if eps == 0
        min_adj_segment_theta = acos(1 - 2/lambda^2);
    else
        min_adj_segment_theta = 0;
    end
    segment_splits = calc_segment_splits(min_adj_segment_theta, lambda, eps, ...
        step_size, min_segment_splits);
    
    % Iterative Process
    for t_n_plus_1 = 2:(n_steps + 1)
        % Get current point and instantiate the pre-calculated range of angles
        % for adjacent segment
        z = points(t_n_plus_1 - 1);
        ranges = [0, 2*pi];
    
        % Determine the directions we can go without violating
        % quasi-geodesic lower bound condition
        for t_i = 1:(t_n_plus_1 - 2) % Note: t_i >= lambda*eps - t_n_plus_1 <=> (t_n_plus_1 - t_i)/lambda - eps < 0
            % Get the i-th previous segment and determine if 
            segment_i = GeodesicSegment(points(t_i), points(t_i + 1));

            % Skip if we are far enough away
            %if segment_i.dist_from_point(z) > (t_n_plus_1 - t_i) / lambda - eps + (2 - 1/lambda)*step_size % thickening + s
                %continue
            %end

            if t_i == t_n_plus_1 - 2
                prev_angle = segment_i.get_angle_with_vertical(step_size);
                range_i = [pi + min_adj_segment_theta + prev_angle, ...
                           3*pi - min_adj_segment_theta + prev_angle];
                ranges = [ranges;range_i];
                continue
            end

            for i = 0:(segment_splits - 1)
                % Determine the s-neighborhood around the segment, which is a
                % neighborhood within which we cannot guarentee stepping into
                % will preserve the quasi-geodesic lower bound
                sub_i = i/segment_splits;
                sub_i_plus_1 = (i + 1)/segment_splits;

                % Calculate proven s. Of course if this s is non-positive,
                % the s-neighborhood will be empty, so we can then skip the
                % following calculations
                lowerBd1 = (1/lambda) * (t_n_plus_1 - 1 - (t_i + sub_i)) * step_size - eps;
                if lowerBd1 >= 0
                    s = acosh((1/lambda - 1)*sinh(step_size)*sinh(lowerBd1) + cosh(step_size + lowerBd1));
                else
                    % If the lower bound is less than 0, then common sense
                    % dictate we need not consider it. However, the issue
                    % lies in the fact that on a given step, s could grow
                    % by up to step_size/lambda + (1 + 1/lambda)step_size.
                    % We thus want to keep our endpoint 
                    s = lowerBd1 + (1 + 1/lambda) * step_size;
                end
                if s <= 0
                    continue
                end
                if t_i == t_n_plus_1 - 3
                    assert(s < step_size, "bruh how... " + s + " v.s. " + step_size + " also " + lowerBd1)
                end

                % We conceptually split the given segment into a number of
                % pieces of equal length
                sub_segment_start = segment_i.travel_from_start(sub_i*step_size);
                sub_segment_end = segment_i.travel_from_start(sub_i_plus_1*step_size);
                sub_segment = GeodesicSegment(sub_segment_start, sub_segment_end);

                % TESTING TO MAKE SURE NEIGHBORHOODS GROW CORRECTLY
                if sub_segment.dist_from_point(z) + step_size <= s
                    "Fatal: Too Close! " + (sub_segment.dist_from_point(z) + step_size) + " vs " + s
                    [t_n_plus_1, t_i, i]
                    continue
                end

                % Getting intersection of neighborhood and step circle
                intersection_i = intersections_of_point_and_segment_ngbhs(z, sub_segment, ...
                step_size, s);

                inside = (sub_segment.dist_from_point(z) <= s);
                range_i = getRangeTn(z, intersection_i, inside, segment_i);
                ranges = [ranges;range_i];
            end
        end
        merged_range = MergeRange(ranges);
        if isempty(merged_range)
            "Stopping at step #" + (t_n_plus_1 - 1)
            points = points(1:t_n_plus_1 - 1);
            break
        end
    
        % Generate the new point in the specified range
        [new_z, phi] = generateTn(z,merged_range,step_size);
        points(t_n_plus_1) = new_z;
        %=&=points

    end
    
    % Plot the path
    figure;
    plot(real(points), imag(points), 'o-');
    xlabel('Real part');
    ylabel('Imaginary part');
    title('Random Walk in the Hyperbolic Upper Half Plane');
    axis equal;

    % Verify quasi-geodesicity
    %=&=points
    %=&=verify_quasigeodesic(points, lambda, eps, step_size, 20)
end

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


% calculate the angle between intersecting points
function rangeTn = getRangeTn(t, intersection, t_inside_boundary, segment)
    % Filter down to 1-2 intersections
    if isempty(intersection)
        rangeTn = [0, 2*pi];
        return
    elseif height(intersection) > 2
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
    %if isempty(mergedRange)
        %double(ranges)
        %error('The ranges do not overlap.');
    %end
end

% generate bounded t_n using randomization
function [tn, phi] = generateTn(t,range,step_size) % ?might be able to generate "geodesic"
    
    % Generate random value within size of range for uniform distribution
    range_size = 0;
    for i = 1:size(range, 1)
        range_size = range_size + (range(i, 2) - range(i, 1));
    end

    % Get phi via random value
    random_value = range_size * rand();
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
    tn = get_point_along_direction(t, phi + pi/2, step_size);
end

function split_val = calc_segment_splits(adj_segment_theta, lambda, eps, ...
    step_size, min_segment_splits)
    if eps == 0
        h = acosh(cosh(step_size)^2 - sinh(step_size)^2*cos(adj_segment_theta));

        lowerBd1 = (1/lambda) * step_size;
        lowerBd2 = 2*lowerBd1; % lowerBd2 = (1/lambda) * 2 * step_size = 2*lowerBd1
        s = acosh((1/lambda - 1)*sinh(step_size)*sinh(lowerBd1) + cosh(step_size + lowerBd1)); 
        thickening = s - lowerBd2;
        max_possible_backstep = lambda * (h + step_size - thickening) - 3*step_size;
        analytical_val = ceil(step_size/max_possible_backstep);
        split_val = max(min_segment_splits, analytical_val);
    else
        split_val = min_segment_splits;
    end
end

% This function validates that the given input parameters
function validate_parameter_conditions(lambda, eps, step_size, n_steps, min_segment_splits)
    if eps < 0
        error("The epsilon value provided (" + eps + ") is negative, but " + ...
            "epsilon must be non-negative.")
    end
    if lambda < 1
        error("The lambda value provided (" + lambda + ") is less than 1, " + ...
            "but lambda must be greater than 1.")
    end
    if eps ~= 0 && step_size > eps
        error("Given epsilon value (" + eps + ") is greater than the give " + ...
            "step size value (" + step_size + "). If epsilon is non-zero, " + ...
            "it is required for accurate results that the provided step " + ...
            "size is smaller than epsilon.")
    end
    if ~is_integer(n_steps)
        error("The provided number of steps to perform (" + n_steps + ") is " + ...
              "not a whole number. Please provide a whole number value greater " + ...
              "than or equal to 2.")
    end
    if n_steps < 2
        error("The provided number of steps to perform (" + n_steps + ") is " + ...
            "less than 2. Please provide a value greater than or equal to 2.")
    end
    if ~is_integer(min_segment_splits)
        error("The provided minimum number of segment splits (" + min_segment_splits + ...
            ") is not a whole number. Please provide a positive integer value.")
    end
    if min_segment_splits < 1
        error("The provided minimum number of segment splits (" + min_segment_splits + ...
            ") is less than one. Please provide a positive integer value.")
    end
end