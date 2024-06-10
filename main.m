
%random_walk_hyperbolic(5, 0.1);

points = get_points_on_circle(1i, 1, 360);
plot(real(points), imag(points))

%x = calc_quasigeodesic_stability_constant(2, 1, 3)