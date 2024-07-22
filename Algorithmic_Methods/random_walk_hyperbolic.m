
% TO-DO: Make z not bad
function [start_points, end_points] = random_walk_hyperbolic(lambda, eps, ...
    step_size, target_length, min_segment_splits, jump_chance, force_beginning_jump)
    % This function performs a random walk in the hyperbolic upper half-plane
    % target_length: The length of the geodesic
    % step_size: Fixed hyperbolic distance to move at each step

    SAFE_RADIUS = 12;
    jump_samplings = min(floor(eps/step_size), 20);

    % Validate that the inputted parameters
    validate_parameter_conditions(lambda, eps, step_size, target_length, min_segment_splits)

    % Initialize the starting point
    z = 1i;  % Start at (0, 0.1) to avoid the real axis (decimal precision gets bad down there)
    
    % Pre-allocate array for storing points of quasi-geodesic
    max_segments = ceil(lambda/step_size * (target_length + eps)) + 1;
    start_points = zeros(max_segments, 1);
    end_points = zeros(max_segments, 1);

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
    for t_n = 1:max_segments
        t_n
        % ====================== EPSILON JUMPING ======================
        % (Determining Start Point of New Segment)
        % Discontinuously jump...
        % =============================================================
        if eps ~= 0 && (jumps(t_n) || (t_n == 1 && force_beginning_jump))
            % 1. Determine if segment can have influence (i.e. if the eps ball
            % pokes into maximum estimate for lower bound or pokes out of minimum
            % estimate for upper bound)
            ti_inside_eps_ball = [];
            for t_i = 1:(t_n - 1)
                segment_i = GeodesicSegment(start_points(t_i), end_points(t_i));
                dist_to_seg = segment_i.dist_from_point(z);
                max_lower_bound = (t_n - t_i)*step_size/lambda - eps; % This can be negative but later logic handles that
                min_upper_bound = lambda*(t_n - (t_i + 1))*step_size + eps;
                if (dist_to_seg - eps < max_lower_bound && max_lower_bound > 0) ...
                   || (dist_to_seg + eps > min_upper_bound)
                    ti_inside_eps_ball = [ti_inside_eps_ball, t_i];
                end
            end

            % 2. Sample jump circles with radius between 0 and... TODO:
            % Finish the comment lol
            if isempty(ti_inside_eps_ball)
                radius = acosh( (cosh(eps) - 1)*rand() + 1); % Inverse CDF w/ PDF as sinh(r)/(cosh(eps) - 1)
                angle = 2*pi*rand();
                new_z = get_point_along_direction(z, angle, radius);
                start_points(t_n) = new_z;

                %display_epsilon_jump_1(z, new_z, eps, t_n, start_points, end_points)
            else
                jump_ranges_in_eps_circle = zeros(0, 4);
                total_weight = 0;

                jump_size_range = linspace(0, eps, jump_samplings + 1);
                for jump_size = jump_size_range(2:jump_samplings + 1) % Skip jump_size = 0
                    ranges_for_jump_circle = [0, 2*pi];
                    jump_circle_unusable = false;
                    for t_i = ti_inside_eps_ball
                        if jump_circle_unusable
                            break
                        end
                        for i = 0:(segment_splits - 1)
                            % Get the percentage along the segment we want to make
                            % our first cut (value needed for following calculations)
                            sub_i = i/segment_splits;
            
                            % Calculate quasi-geodesic lower bound and
                            % upper bound
                            dist_in_R = (t_n - (t_i + sub_i)) * step_size;
                            lower_bound = (1/lambda) * dist_in_R - eps;
                            upper_bound = lambda * dist_in_R + eps;
            
                            % Retrieve cahsed data for subsegment
                            sub_segment_start = sub_segment_points(t_i, i + 1);
                            sub_segment_end = sub_segment_points(t_i, i + 2);
                            sub_segment = GeodesicSegment(sub_segment_start, sub_segment_end);
                            e1 = sub_segment_points_transformed(t_i, 2*i + 1);
                            e2 = sub_segment_points_transformed(t_i, 2*i + 2);
                            a = sub_segment_points_abcd_values(t_i, 4*i + 1);
                            b = sub_segment_points_abcd_values(t_i, 4*i + 2);
                            c = sub_segment_points_abcd_values(t_i, 4*i + 3);
                            d = sub_segment_points_abcd_values(t_i, 4*i + 4);

                            dist_to_subseg = sub_segment.dist_from_point(z);

                            % I. Lower Bound
                            if lower_bound >= 0
                                % If jump circle is complete engulfed by
                                % lower bound, the entire circle is not
                                % safe to jump to
                                if dist_to_subseg + jump_size <= lower_bound
                                    jump_circle_unusable = true;
                                    break
                                end

                                lb_intersections = intersections_of_point_and_segment_ngbhs(z, ...
                                e1, e2, jump_size, lower_bound, a, b, c, d);
                                
                                if ~isempty(lb_intersections)
                                    range_i = getRangeTn(z, lb_intersections, sub_segment, ...
                                                         jump_size, lower_bound, true);
                                    ranges_for_jump_circle = [ranges_for_jump_circle; range_i];
                                end
                            end
                            
                            % II. Upper Bound
                            % If jump circle is completely outside of the 
                            % upper bound, the entire circle is not safe to 
                            % jump to
                            if jump_size >= min_upper_bound + upper_bound
                                jump_circle_unusable = true;
                                break
                            end
            
                            % Get the intersection of lower/upper bound 
                            % neighborhood and jump circle
                            if sub_segment.dist_from_point(z) + jump_size >= upper_bound
                                ub_intersections = intersections_of_point_and_segment_ngbhs(z, ...
                                e1, e2, jump_size, upper_bound, a, b, c, d);
                
                                % Get range of angles to not cross into lower/upper 
                                % bound
                                if ~isempty(ub_intersections)
                                    range_i = getRangeTn(z, ub_intersections, sub_segment, ...
                                                         jump_size, upper_bound, false);
                                    ranges_for_jump_circle = [ranges_for_jump_circle; range_i];
                                end
                            end
                        end
                    end
                    if ~jump_circle_unusable
                        merged_ranges = merge_ranges(ranges_for_jump_circle);
                        for i = 1:height(merged_ranges)
                            theta_1 = merged_ranges(i, 1);
                            theta_2 = merged_ranges(i, 2);
                            weight = (theta_2 - theta_1) * 2 * pi * sinh(jump_size);
                            jump_ranges_in_eps_circle = [jump_ranges_in_eps_circle; 
                                                         theta_1, theta_2, jump_size, weight];
                            total_weight = total_weight + weight;
                        end
                    end
                end
                % JUMP!
                if isempty(jump_ranges_in_eps_circle)
                    start_points = start_points(1:t_n-1);
                    end_points = end_points(1:t_n-1);
        
                    disp("Nowhere safe to jump. Stopping at step #" + t_n)
                    break  
                end
                new_z = generateJumpPoint(z, jump_ranges_in_eps_circle, ...
                                                     total_weight);
                start_points(t_n) = new_z;

                % display_epsilon_jump_2(z, new_z, lambda, eps, step_size, t_n, segment_splits, ti_inside_eps_ball, ...
                %    jump_ranges_in_eps_circle, jump_size_range, start_points, end_points, ...
                %     sub_segment_points, sub_segment_points_transformed, sub_segment_points_abcd_values)
            end
        else
            start_points(t_n) = z;
        end

        % ===================== CHOOSING NEXT POINT =====================
        % (Determining End Point of New Segment)
        % Determine the directions we can travel along our step circle 
        % without violating quasi-geodesic lower bound condition
        % ===============================================================
        z = start_points(t_n);

        % Instaniate the pre-calculated range of angles for adjacent segment
        allowed_ranges = [0, 2*pi];
        % Get angle range we can step along
        for t_i = 1:(t_n - 1)
            % Get the i-th previous segment
            %segment_i = GeodesicSegment(start_points(t_i), end_points(t_i));

            % Skip if we are far enough away
            %if segment_i.dist_from_point(z) > (t_n_plus_1 - t_i) / lambda - eps + (2 - 1/lambda)*step_size % thickening + s
                %continue
            %end

            % SPECIAL CASE: If the segment we are looking at is the last
            % segment in the path. In this case, since step_size <= eps, this
            % segment does not produce a lower bound. However, if eps = 0,
            % it will still produce a lower bound, so we use the
            % pre-calculated angle range above
            if t_i == t_n - 1
                if eps == 0
                    % Shift pre-calculated range to be within the
                    % directional frame along the segment of our 
                    % current point z
                    segment_i = GeodesicSegment(start_points(t_i), end_points(t_i));
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
                    if s <= 1e-10
                        continue
                    end
                end
               
                % Retrieve cahsed data
                sub_segment_start = sub_segment_points(t_i, i + 1);
                sub_segment_end = sub_segment_points(t_i, i + 2);
                sub_segment = GeodesicSegment(sub_segment_start, sub_segment_end);
                e1 = sub_segment_points_transformed(t_i, 2*i + 1);
                e2 = sub_segment_points_transformed(t_i, 2*i + 2);
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
                    range_i = getRangeTn(z, intersection_i, sub_segment, ...
                                         step_size, s, true);
                    allowed_ranges = [allowed_ranges; range_i];
                end
            end
        end
        % Merge all allowed ranges by taking their intersection, and pick
        % a direction from this range to step along
        merged_range = merge_ranges(allowed_ranges);
        if isempty(merged_range)
            start_points = start_points(1:t_n-1);
            end_points = end_points(1:t_n-1);

            allowed_ranges
            disp("No overlapping step ranges. Stopping at step #" + t_n)
            break
        end 
        
        % Generate the new point
        [new_z, phi] = generateTn(z,merged_range,step_size);

        %display_random_walk(z, new_z, lambda, eps, step_size, segment_splits, t_n, ...
        %   merged_range, sub_segment_points, sub_segment_points_transformed, ...
        %   sub_segment_points_abcd_values);


        end_points(t_n) = new_z;
        z = new_z;

        % Add new point and update cached/stored data for optimization
        [sub_segment_points, ...
         sub_segment_points_transformed, ...
         sub_segment_points_abcd_values] = ...
         store_new_segment_data(t_n, GeodesicSegment(start_points(t_n), end_points(t_n)), ...
                              step_size, segment_splits, ...
                              sub_segment_points, ...
                              sub_segment_points_transformed, ...
                              sub_segment_points_abcd_values);

        % Normalize into safe zone if necessary
        if dist_H(start_points(t_n), 1i) >= SAFE_RADIUS || dist_H(end_points(t_n), 1i) >= SAFE_RADIUS
            [max_length, p1, p2] = find_maximizing_points(start_points, end_points, t_n);
            if max_length >= 2*SAFE_RADIUS
                "Breaking: Safe Limit Reached"
                break
            end
            p1_p2_midpoint = GeodesicSegment(p1, p2).get_midpoint();
            % Calculate a and b for FLT (we always have c=0 and d=1)
            a = 1/imag(p1_p2_midpoint);
            b = 0 - a*real(p1_p2_midpoint);
            % Transform all points and then recalculate cached data
            for i = 1:t_n
                start_points(i) = a*start_points(i) + b;
                end_points(i) = a*end_points(i) + b;
            end
            [sub_segment_points, ...
             sub_segment_points_transformed, ...
             sub_segment_points_abcd_values] = recalc_stored_data(start_points, end_points, ...
                                                                  t_n, step_size, segment_splits, ...
                                                                  sub_segment_points, ...
                                                                  sub_segment_points_transformed, ...
                                                                  sub_segment_points_abcd_values);
        end
    end

    % Cut off quasi-geodesic at target length
    num_points = length(end_points);
    last_index = -1;
    % Loop backwards to find last point within eps of the target length
    % from the start
    if eps ~= 0
        for i = num_points:-1:1
            dist = dist_H(end_points(i), start_points(1));
            if abs(dist - target_length) <= eps
                last_index = i;
                break
            end
        end
    else
        prev_dist = -1;
        for i = (num_points - 1):-1:1
            p_matrix = [end_points(i), start_points(i + 1);
                        start_points(i), end_points(i)];
            for j = 1:2
                p1 = p_matrix(j, 1);
                p2 = p_matrix(j, 2); 
                dist = dist_H(start_points(1), p1);
                if dist <= target_length && prev_dist > target_length
                    prev_dist = dist;
                    g = GeodesicSegment(p1, p2);
                    u = real(start_points(1));
                    v = imag(start_points(1));
                    c = g.get_center_on_real_line();
                    r = g.get_radius_from_center();
                    y_lb = min(imag(p1), imag(p2));
                    y_ub = max(imag(p1), imag(p2));
                    last_point = find_intersection_of_circles( ...
                                         u, v*cosh(target_length), v*sinh(target_length), ...
                                         c, 0, r, ...
                                         y_lb, y_ub);
                    if ~isempty(last_point)
                        if j == 2
                            start_points(i) = last_point;
                        end
                        end_points(i) = last_point;
                        last_index = i;
                        break
                    end
                end
                prev_dist = dist;
            end
        end
    end
    assert(last_index ~= -1, "ERROR: No point within epsilon of the target " + ...
        "length was found!")

    % Update start points and end points to the desired cutoff index
    start_points = start_points(1:last_index);
    end_points = end_points(1:last_index);
    
    % Plot the path
    plot_quasigeodesic(start_points, end_points)
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

function [max_length, p1, p2] = find_maximizing_points(start_points, end_points, t_n)
    max_length = 0;
    for i = 1:t_n
        for j = 1:t_n
            p_matrix = [
                start_points(i), start_points(j);
                start_points(i), end_points(j);
                end_points(i), start_points(j);
                end_points(i), end_points(j);
            ];
            for k = 1:4 
                candidate_p1 = p_matrix(k, 1);
                candidate_p2 = p_matrix(k, 2);
                dist = dist_H(candidate_p1, candidate_p2);
                if dist > max_length
                    max_length = dist;
                    p1 = candidate_p1;
                    p2 = candidate_p2;
                end
            end
        end
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