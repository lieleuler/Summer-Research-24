
function display_epsilon_jump_1(z, new_z, eps, t_n, start_points, end_points)
    hold on

    % Segments + Points
    for t_i = 1:(t_n - 1)
        segment_i = GeodesicSegment(start_points(t_i), end_points(t_i));
        segment_i.plot(100, "k")
        plot(real(start_points(t_i)), imag(start_points(t_i)), "o", "MarkerSize", 4, "MarkerEdgeColor", "k", "MarkerFaceColor", "k");
    end

    % Circle
    plot_circular_arc(z, eps, 0, 2*pi, "m")

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