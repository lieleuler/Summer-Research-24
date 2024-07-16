
function plot_circular_arc(center, radius, theta_1, theta_2, color)
    thetas = linspace(theta_1, theta_2, 100);
    points = zeros(1, 100);
    for i = 1:100
        theta = thetas(i);
        points(i) = get_point_along_direction(center, theta, radius);
    end
    plot(real(points), imag(points), color)
end