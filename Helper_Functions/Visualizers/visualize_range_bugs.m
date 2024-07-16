
function visualize_range_bug(z, t_n, lambda, eps, step_size, segment_splits, min_adj_segment_theta, ...
    start_points, end_points, sub_segment_points, sub_segment_points_transformed, sub_segment_points_abcd_values)

    % ===================== CHOOSING NEXT POINT =====================
    % (Determining End Point of New Segment)
    % Determine the directions we can travel along our step circle 
    % without violating quasi-geodesic lower bound condition
    % ===============================================================
    z = start_points(t_n);

    % Instaniate the pre-calculated range of angles for adjacent segment
    allowed_ranges = [0, 2*pi];
    corresponding_segs = GeodesicSegment(0);
    corresponding_s = zeros(0);
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
                corresponding_segs(i) = sub_segment;
                corresponding_s(i) = s;
                range_i = getRangeTn(z, intersection_i, sub_segment, ...
                                     step_size, s, true);
                allowed_ranges = [allowed_ranges; range_i];
            end
        end
    end
    
    N = length(corresponding_segs);
    gradient = generate_rainbow_gradient(N);
    for i = 1:N
        curr_seg = corresponding_segs(i);
        curr_s = corresponding_s(i);
        curr_range = allowed_ranges(i);
        
        curr_seg.plot_neighborhood(curr_s, gradient(i));
        plot_circular_arc(z, step_size/2 + step_size/2*i/N, curr_range(1), curr_range(2), gradient(i))
    end
end