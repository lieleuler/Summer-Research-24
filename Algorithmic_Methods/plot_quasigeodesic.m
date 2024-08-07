
function plot_quasigeodesic(start_points, end_points)
    N = length(start_points);
    gradient = generate_rainbow_gradient(N);
    figure;
    hold on
    for i = 1:N
        plot([real(double(start_points(i))), real(double(end_points(i)))], ...
             [imag(double(start_points(i))), imag(double(end_points(i)))], ...
             'o-', "Color", gradient(i, :));
    end
    xlabel('Real part');
    ylabel('Imaginary part');
    title('Random Walk in the Hyperbolic Upper Half Plane');
    axis equal;
end