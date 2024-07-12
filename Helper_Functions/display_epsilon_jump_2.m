
function display_epsilon_jump_2(z, new_z, lambda, eps, step_size, t_n, segment_splits, ti_inside_eps_ball, ...
                jump_ranges_in_eps_circle, jump_size_range, start_points, end_points, ...
                sub_segment_points, sub_segment_points_transformed, sub_segment_points_abcd_values)
    
    hold on
    % Segments
    for t_i = 1:(t_n - 1)
        segment_i = GeodesicSegment(start_points(t_i), end_points(t_i));
        segment_i.plot(100, "k")
        plot(real(start_points(t_i)), imag(start_points(t_i)), "o", "MarkerSize", 4, "MarkerEdgeColor", "k", "MarkerFaceColor", "k");
    end

    % Relevant Neighborhoods
    for t_i = ti_inside_eps_ball
        for i = 0:(segment_splits - 1)
            sub_i = i/segment_splits;

            dist_in_R = (t_n - (t_i + sub_i)) * step_size;
            lower_bound = (1/lambda) * dist_in_R - eps;
            upper_bound = lambda * dist_in_R + eps;
    
            sub_segment_start = sub_segment_points(t_i, i + 1);
            sub_segment_end = sub_segment_points(t_i, i + 2);
            sub_segment = GeodesicSegment(sub_segment_start, sub_segment_end);
            e1 = sub_segment_points_transformed(t_i, 2*(i + 1) - 1);
            e2 = sub_segment_points_transformed(t_i, 2*(i + 1) - 1);
            a = sub_segment_points_abcd_values(t_i, 4*i + 1);
            b = sub_segment_points_abcd_values(t_i, 4*i + 2);
            c = sub_segment_points_abcd_values(t_i, 4*i + 3);
            d = sub_segment_points_abcd_values(t_i, 4*i + 4);

            if lower_bound > 0
                sub_segment.plot_neighborhood(lower_bound, "r");
            end
            if upper_bound > 0
                sub_segment.plot_neighborhood(upper_bound, "b");
            end
        end
    end

    % Jump Circles
    for jump_size = jump_size_range
       plot_circular_arc(z, jump_size, 0, 2*pi, "m")
    end
    
    % Good Ranges
    h = jump_ranges_in_eps_circle
    for i = 1:height(jump_ranges_in_eps_circle)
       theta_1 = jump_ranges_in_eps_circle(i, 1);
       theta_2 = jump_ranges_in_eps_circle(i, 2);
       r = jump_ranges_in_eps_circle(i, 3);
       plot_circular_arc(z, r, theta_1, theta_2, "g")
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