
function intersections = intersections_of_circles(center1, radius1, center2, radius2)
    c1_x = real(center1);
    c1_y = imag(center1);
    c2_x = real(center2);
    c2_y = imag(center2);

    syms x y

    circle1_eq = (x - c1_x)^2 + (y - c1_y)^2 - 2*c1_y*y*(cosh(radius1) - 1) == 0;
    circle2_eq = (x - c2_x)^2 + (y - c2_y)^2 - 2*c2_y*y*(cosh(radius2) - 1) == 0;

    sols = solve(circle1_eq, circle2_eq, [x, y]);

    intersections = [sols.x + i*sols.y];