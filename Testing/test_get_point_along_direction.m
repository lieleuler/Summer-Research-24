z = 2.3 + 0.8i;
angles = linspace(0, 2*pi, 301);

points_x = zeros(length(angles));
points_y = zeros(length(angles));
tic
for i = 1:length(angles)
    circle_point = get_point_along_direction(z, angles(i), 1);
    points_x(i) = real(circle_point);
    points_y(i) = imag(circle_point);
end
toc

plot(points_x, points_y)