
function intersections = intersections_of_point_and_segment_ngbhs(point, ...
    geodesic, p_radius, g_radius)

    [a, b, c, d] = geodesic.find_flt_to_imag_axis();
    [e1, e2] = geodesic.fractional_linear_transform(a, b, c, d).get_endpoints();
    if imag(e1) < imag(e2)
        a1 = imag(e1);
        a2 = imag(e2);
    else
        a1 = imag(e2);
        a2 = imag(e1);
    end

    trans_point = (a*point + b) / (c*point + d);
    u = real(trans_point);
    v = imag(trans_point);

    syms x y
    point_ngbh_eq = (x - u)^2 + (y - v)^2 - 2*v*y*(cosh(p_radius) - 1) == 0;

    % Eq 1: Lower Circle
    assume(y < a1*sech(g_radius))
    bottom_circle_eq = x^2 + (y - a1*cosh(g_radius))^2 == (a1*sinh(g_radius))^2;
    sols1 = solve(point_ngbh_eq, bottom_circle_eq, [x, y]);

    % Eq 2 + 3: Outer Lines
    assume(y, "clear")
    assume(y >= a1*sech(g_radius))
    assumeAlso(y <= a2*sech(g_radius))
    right_line = y == csch(g_radius)*(x - a1*tanh(g_radius)) + a1*sech(g_radius);
    left_line = y == csch(g_radius)*(-x - a1*tanh(g_radius)) + a1*sech(g_radius);
    sols2 = solve(point_ngbh_eq, right_line, [x, y]);
    sols3 = solve(point_ngbh_eq, left_line, [x, y]);

    % Eq 4: Top Circle
    assume(y, "clear")
    assume(y > a2*sech(g_radius))
    top_circle_eq = x^2 + (y - a2*cosh(g_radius))^2 == (a2*sinh(g_radius))^2;
    sols4 = solve(point_ngbh_eq, top_circle_eq, [x, y]);

    intersections = [sols1.x + 1i*sols1.y;
                     sols2.x + 1i*sols2.y;
                     sols3.x + 1i*sols3.y;
                     sols4.x + 1i*sols4.y];
    
end