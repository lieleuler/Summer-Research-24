
function points = random_walk_hyperbolic(lambda, eps, ...
    step_size, target_length, min_segment_splits, jump_chance)
    % This function performs a random walk in the hyperbolic upper half-plane
    % target_length: The length of the geodesic
    % step_size: Fixed hyperbolic distance to move at each step

    % Validate that the inputted parameters
    validate_parameter_conditions(lambda, eps, step_size, target_length, min_segment_splits)

    % Initialize the starting point
    z = 1i;  % Start at (0, 1) to avoid the real axis
    
    % Pre-allocate array for storing points of quasi-geodesic
    max_segments = ceil(lambda/step_size * (target_length + eps));
    points = zeros(max_segments + 1, 1);
    points(1) = z;

    % Pre-calculate the angle range allowed by the lower bound of the
    % segment adjacent to the point (this only matters if eps = 0)
    if eps == 0
        min_adj_segment_theta = acos(1 - 2/lambda^2);
    else
        min_adj_segment_theta = 0;
    end

    % Pre-calculate the amount of pieces to split each segment into
    % during step_circle-segment_neighborhood intersection calculation
    segment_splits = calc_segment_splits(min_adj_segment_theta, lambda, eps, ...
        step_size, min_segment_splits);

    % Initialize pre-allocted caching tables (these appear to greatly
    % improve runtime)
    sub_segment_points = zeros(max_segments + 1, segment_splits + 1);
    sub_segment_points_transformed = zeros(max_segments + 1, 2*segment_splits);
    sub_segment_points_abcd_values = zeros(max_segments + 1, 4*segment_splits);

    % Initialize jumping array
    random_numbers = rand(max_segments, 1);
    jumps = random_numbers <= jump_chance;
    
    % Iterative Process
    t_n_plus_1 = 1;
    while true
        t_n_plus_1 = t_n_plus_1 + 1;
        t_n = t_n_plus_1 - 1;
        % Instaniate the pre-calculated range of angles for adjacent segment
        rows = segment_splits*(t_n_plus_1 - 2) + 1;
        allowed_ranges = [0, 2*pi];
    
        % Determine the directions we can travel along our step circle 
        % without violating quasi-geodesic lower bound condition
        for t_i = 1:(t_n_plus_1 - 2) % Note: t_i >= lambda*eps - t_n_plus_1 <=> (t_n_plus_1 - t_i)/lambda - eps < 0
            % Get the i-th previous segment
            segment_i = GeodesicSegment(points(t_i), points(t_i + 1));

            % Skip if we are far enough away
            %if segment_i.dist_from_point(z) > (t_n_plus_1 - t_i) / lambda - eps + (2 - 1/lambda)*step_size % thickening + s
                %continue
            %end

            % SPECIAL CASE: If the segment we are looking at is the last
            % segment in the path. In this case, since step_size <= eps, this
            % segment does not produce a lower bound. However, if eps = 0,
            % it will still produce a lower bound, so we use the
            % pre-calculated angle range above
            if t_i == t_n_plus_1 - 2
                if eps == 0
                    % Shift pre-calculated range to be within the
                    % directional frame along the segment of our 
                    % current point z
                    prev_angle = segment_i.get_angle_with_vertical(step_size);
                    range_i = [pi + min_adj_segment_theta + prev_angle, ...
                               3*pi - min_adj_segment_theta + prev_angle];
                    allowed_ranges = [allowed_ranges; range_i];
                end
                continue
            end

            % Split segment into pieces of equal length, and determine 
            % the s-neighborhood around the segment, which is a
            % neighborhood within which we cannot guarentee stepping 
            % into will preserve the quasi-geodesic lower bound
            for i = 0:(segment_splits - 1)
                % Get the percentage along the segment we want to make
                % our first cut (value needed for following calculations)
                sub_i = i/segment_splits;

                % Calculate proven s. Of course if this s is non-positive,
                % the s-neighborhood will be empty, so we can then skip the
                % following calculations if this is the case
                lowerBd1 = (1/lambda) * (t_n - (t_i + sub_i)) * step_size - eps;
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

                % Retrieve cahsed data
                sub_segment_start = sub_segment_points(t_i, i + 1);
                sub_segment_end = sub_segment_points(t_i, i + 2);
                sub_segment = GeodesicSegment(sub_segment_start, sub_segment_end);
                e1 = sub_segment_points_transformed(t_i, 2*(i + 1) - 1);
                e2 = sub_segment_points_transformed(t_i, 2*(i + 1) - 1);
                a = sub_segment_points_abcd_values(t_i, 4*i + 1);
                b = sub_segment_points_abcd_values(t_i, 4*i + 2);
                c = sub_segment_points_abcd_values(t_i, 4*i + 3);
                d = sub_segment_points_abcd_values(t_i, 4*i + 4);


                %TESTING TO MAKE SURE NEIGHBORHOODS GROW CORRECTLY
                %if sub_segment.dist_from_point(z) + step_size <= s
                    %"Fatal: Too Close! " + (sub_segment.dist_from_point(z) + step_size) + " vs " + s
                    %[t_n_plus_1, t_i, i]
                    %continue
                %end

                % Get the intersection of neighborhood and step circle
                intersection_i = intersections_of_point_and_segment_ngbhs(z, ...
                e1, e2, step_size, s, a, b, c, d);

                % If there is an intersection, we calculate the range of
                % directions - represented by angles - we can travel along
                % the step circle without entering the lower bound
                if ~isempty(intersection_i)
                    %inside = (sub_segment.dist_from_point(z) <= s);
                    range_i = getRangeTn(z, intersection_i, sub_segment, ...
                                         step_size, s);
                    allowed_ranges = [allowed_ranges; range_i];
                end
            end
        end
        % Merge all allowed ranges by taking their intersection, and pick
        % a direction from this range to step along
        merged_range = merge_ranges(allowed_ranges);
        if isempty(merged_range)
            allowed_ranges
            "Stopping at step #" + (t_n)
            points = points(1:t_n);
            break
        end

        % Generate the new point
        [new_z, phi] = generateTn(z,merged_range,step_size);
        length_end_to_end = dist_H(points(1), new_z);
        if length_end_to_end >= target_length
            new_z = get_point_along_direction(z, phi, length_end_to_end - target_length);
            points(t_n_plus_1) = new_z;
            points = points(1:t_n);
            break
        end
        points(t_n_plus_1) = new_z;
        z = new_z;


        % Add new point and update cached/stored data for optimization
        [sub_segment_points_transformed, ...
         sub_segment_points_abcd_values] = ...
        store_new_segment_data(t_n, GeodesicSegment(points(t_n), new_z), ...
                              step_size, segment_splits, ...
                              sub_segment_points_transformed, ...
                              sub_segment_points_abcd_values);

    end
    
    % Plot the path
    figure;
    plot(real(points), imag(points), 'o-');
    xlabel('Real part');
    ylabel('Imaginary part');
    title('Random Walk in the Hyperbolic Upper Half Plane');
    axis equal;
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
function validate_parameter_conditions(lambda, eps, step_size, target_length, min_segment_splits)
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
    if target_length <= 0
         error("The provided target length (" + target_length + ") is " + ...
              "non-positive. Please provide a positive whole number value")
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