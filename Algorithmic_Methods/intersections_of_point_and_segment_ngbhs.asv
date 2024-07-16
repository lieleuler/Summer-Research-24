function intersections = intersections_of_point_and_segment_ngbhs(point, ...
    e1, e2, p_radius, g_radius, a, b, c, d)

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

    q = 2*v*(cosh(p_radius) - 1);

    intersections = zeros(0, 1);

    % disp([a1, a2])
    % disp(u)
    % disp(v)

    % Eq 1: Lower Circle
    h = a1*cosh(g_radius);
    R = a1*sinh(g_radius);
    A = 2*h - 2*v - q;
    B = R^2 - h^2 + v^2 + u^2;
    a_0 = A^2 + 4*u^2;
    b_0 = 2*(A*B - 4*u^2*h);
    c_0 = B^2 + 4*u^2*(h^2 - R^2);
    discriminant = b_0^2 - 4*a_0*c_0;
    if discriminant >= 0
        y_1 = (-b_0 + sqrt(discriminant)) / (2*a_0);
        y_2 = (-b_0 - sqrt(discriminant)) / (2*a_0);
        if y_1 < a1*sech(g_radius)
            x_1 = sqrt(R^2 - y_1^2 + 2*h*y_1 - h^2);
            if abs(x_1) < a1*tanh(g_radius)
                if abs((x_1 - u)^2 + (y_1 - v)^2 - q*y_1) > abs((-x_1 - u)^2 + (y_1 - v)^2 - q*y_1)
                    x_1 = -x_1;
                end
                sol = (d*(x_1 + 1i*y_1) - b)/(-c*(x_1 + 1i*y_1) + a);
                intersections = [intersections; sol];
            end
        end
        if y_2 < a1*sech(g_radius)
            x_2 = sqrt(R^2 - y_2^2 + 2*h*y_2 - h^2);
            if abs(x_2) < a1*tanh(g_radius)
                if abs((x_2 - u)^2 + (y_2 - v)^2 - q*y_2) > abs((-x_2 - u)^2 + (y_2 - v)^2 - q*y_2)
                    x_2 = -x_2;
                end 
                sol = (d*(x_2 + 1i*y_2) - b)/(-c*(x_2 + 1i*y_2) + a);
                intersections = [intersections; sol];
            end
        end
    end

    % Eq 2: Upper Circle
    h = a2*cosh(g_radius);
    R = a2*sinh(g_radius);
    A = 2*h - 2*v - q;
    B = R^2 - h^2 + v^2 + u^2;
    a_0 = A^2 + 4*u^2;
    b_0 = 2*(A*B - 4*u^2*h);
    c_0 = B^2 + 4*u^2*(h^2 - R^2);
    discriminant = b_0^2 - 4*a_0*c_0;
    if discriminant >= 0
        y_1 = (-b_0 + sqrt(discriminant)) / (2*a_0);
        y_2 = (-b_0 - sqrt(discriminant)) / (2*a_0);
        if y_1 > a2*sech(g_radius)
            x_1 = sqrt(R^2 - y_1^2 + 2*h*y_1 - h^2);
            if abs(x_1) > a2*tanh(g_radius) || y_1 >= a2*(2*cosh(g_radius) - sech(g_radius))
                if abs((x_1 - u)^2 + (y_1 - v)^2 - q*y_1) > abs((-x_1 - u)^2 + (y_1 - v)^2 - q*y_1)
                    x_1 = -x_1;
                end
                sol = (d*(x_1 + 1i*y_1) - b)/(-c*(x_1 + 1i*y_1) + a);
                intersections = [intersections; sol];
            end
        end
        if y_2 > a2*sech(g_radius)
            x_2 = sqrt(R^2 - y_2^2 + 2*h*y_2 - h^2);
            if abs(x_2) > a2*tanh(g_radius) || y_2 >= a2*(2*cosh(g_radius) - sech(g_radius))
                if abs((x_2 - u)^2 + (y_2 - v)^2 - q*y_2) > abs((-x_2 - u)^2 + (y_2 - v)^2 - q*y_2)
                    x_2 = -x_2;
                end
                sol = (d*(x_2 + 1i*y_2) - b)/(-c*(x_2 + 1i*y_2) + a);
                intersections = [intersections; sol];
            end
        end
    end


    % Eq 3: Right Outer Line
    m = csch(g_radius);
    a_0 = 1/m^2 + 1;
    b_0 = -2*u/m - 2*v - q;
    c_0 = u^2 + v^2;

    discriminant = b_0^2 - 4*a_0*c_0;
    if discriminant >= 0
        y_1 = (-b_0 + sqrt(discriminant)) / (2*a_0);
        y_2 = (-b_0 - sqrt(discriminant)) / (2*a_0);
        if a1*sech(g_radius) <= y_1 && y_1 <= a2*sech(g_radius)
            x_1 = y_1/m;
            if a1*tanh(g_radius) <= abs(x_1) && abs(x_1) <= a2*tanh(g_radius)
                sol = (d*(x_1 + 1i*y_1) - b)/(-c*(x_1 + 1i*y_1) + a);
                intersections = [intersections; sol];
            end
        end
        if a1*sech(g_radius) <= y_2 && y_2 <= a2*sech(g_radius)
            x_2 = y_2/m;
            if a1*tanh(g_radius) <= abs(x_2) && abs(x_2) <= a2*tanh(g_radius)
                sol = (d*(x_2 + 1i*y_2) - b)/(-c*(x_2 + 1i*y_2) + a);
                intersections = [intersections; sol];
            end
        end
    end

    % Eq 4: Left Outer Line
    m = -csch(g_radius);
    a_0 = 1/m^2 + 1;
    b_0 = -2*u/m - 2*v - q;
    c_0 = u^2 + v^2;

    discriminant = b_0^2 - 4*a_0*c_0;
    if discriminant >= 0
        y_1 = (-b_0 + sqrt(discriminant)) / (2*a_0);
        y_2 = (-b_0 - sqrt(discriminant)) / (2*a_0);
        if a1*sech(g_radius) <= y_1 && y_1 <= a2*sech(g_radius)
            x_1 = y_1/m;
            if a1*tanh(g_radius) <= abs(x_1) && abs(x_1) <= a2*tanh(g_radius)
                sol = (d*(x_1 + 1i*y_1) - b)/(-c*(x_1 + 1i*y_1) + a);
                intersections = [intersections; sol];
            end
        end
        if a1*sech(g_radius) <= y_2 && y_2 <= a2*sech(g_radius)
            x_2 = y_2/m;
            if a1*tanh(g_radius) <= abs(x_2) && abs(x_2) <= a2*tanh(g_radius)
                sol = (d*(x_2 + 1i*y_2) - b)/(-c*(x_2 + 1i*y_2) + a);
                intersections = [intersections; sol];
            end
        end
    end

    % Transform Intersections Back
    %intersections = (d*intersections - b)./(-c*intersections + a);
end