
function int = find_intersection_of_circles(x1, y1, r1, x2, y2, r2, y_lb, y_ub)
    A = (r1^2 - r2^2 - x1^2 - y1^2 + x2^2 + y2^2)/(2*(x2 - x1)) - x1;
    B = (y1 - y2)/(x2 - x1);

    a = (B^2 + 1);
    b = 2*(A*B - y1);
    c = A^2 + y1^2 - r1^2;
    
    discriminant = b^2 - 4*a*c;
    if discriminant < 0
        "m"
        int = [];
        return
    end
    y_1 = (-b + sqrt(discriminant))/(2*a);
    y_2 = (-b - sqrt(discriminant))/(2*a);
    if y_lb <= y_1 && y_1 <= y_ub
        y = y_1;
    elseif y_lb <= y_2 && y_2 <= y_ub
        y = y_2;
    else
        "n"
        int = [];
        return
    end

    x_candidate_1 = x1 + sqrt(r1^2 - (y - y1)^2);
    x_candidate_2 = x1 - sqrt(r1^2 - (y - y1)^2);
    candidate_1_error = abs((x_candidate_1 - x2)^2 + (y - y2)^2 - r2^2)
    candidate_2_error = abs((x_candidate_2 - x2)^2 + (y - y2)^2 - r2^2)
    if candidate_1_error <= candidate_2_error
        int = x_candidate_1 + y*1i;
    else
        int = x_candidate_2 + y*1i;
    end
end

