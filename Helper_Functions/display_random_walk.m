
function display_random_walk(z, new_z, lambda, eps, step_size, segment_splits, t_n, ...
    merged_range, sub_segment_points, sub_segment_points_transformed, sub_segment_points_abcd_values)
    hold on

    % Segments
    for t_i = 1:(t_n - 1)
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
            s = max(s, 0);

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

            if s > 0
                sub_segment.plot(100, "k")
            else
                sub_segment.plot(100, "m")
            end

            if t_i ~= t_n - 1
                sub_segment.plot_neighborhood(s, "r")
            end
        end
    end

    % Step Circle
    plot_circular_arc(z, step_size, 0, 2*pi, "b")
    for i = 1:height(merged_range)
        plot_circular_arc(z, step_size, merged_range(i, 1), ...
                          merged_range(i, 2), "g")
    end

    % New Point
    pink = [1, 0.4, 0.8];
    plot(real(new_z), imag(new_z), "p", "MarkerSize", 7, "MarkerEdgeColor", pink, "MarkerFaceColor", pink);

    % Pausing
    while true
        key = waitforbuttonpress;
        if key == 1
            break
        end
    end

    % Clear
    clc
    clf
end