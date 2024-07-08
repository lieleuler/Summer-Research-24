
% This function finds the intersections of a segment s with the bisector
% of two segments s_j and s_k. It first calculates all 14 possible
% intersections between s and the bisector, and then rules out candidates
% based on their coordinates
function points = find_bisector_intersections(s, s_j, s_k)

    % Step 1: Get candidate points from point-point bisectors. 
    points = zeros(0, 1);

    [ej_1, ej_2] = s_j.get_endpoints();
    [ek_1, ek_2] = s_k.get_endpoints();

    new_points = point_point_bisector_intersections(s, ej_1, ek_1, s_j, s_k, 1, 1);
    points = [points; new_points];

    new_points = point_point_bisector_intersections(s, ej_1, ek_2, s_j, s_k, 1, 3);
    points = [points; new_points];

    new_points = point_point_bisector_intersections(s, ej_2, ek_1, s_j, s_k, 3, 1);
    points = [points; new_points];

    new_points = point_point_bisector_intersections(s, ej_2, ek_2, s_j, s_k, 3, 3);
    points = [points; new_points];

    % Step 2: Transform s_j to imaginary axis
    [a, b, c, d] = s_j.find_flt_to_imag_axis();
    trans_sj_1 = s_j.fractional_linear_transform(a, b, c, d);
    trans_sk_1 = s_k.fractional_linear_transform(a, b, c, d);
    trans_s_1 = s.fractional_linear_transform(a, b, c, d);
    [trans_ek_1, trans_ek_2] = trans_sk_1.get_endpoints();

    trans_points_1 = zeros(0, 1);

    new_points = line_line_bisectors_intersection(trans_s_1, trans_sj_1, trans_sk_1);
    trans_points_1 = [trans_points_1; new_points];

    new_points = point_line_bisector_intersection(trans_s_1, trans_ek_1, trans_sj_1, trans_sk_1, 1);
    trans_points_1 = [trans_points_1; new_points];

    new_points = point_line_bisector_intersection(trans_s_1, trans_ek_2, trans_sj_1, trans_sk_1, 3);
    trans_points_1 = [trans_points_1; new_points];

    points_1 = (d*trans_points_1 - b)./(-c*trans_points_1 + a);

    % Step 3: Transform s_k to imaginary axis
    [a, b, c, d] = s_k.find_flt_to_imag_axis();
    trans_sj_2 = s_j.fractional_linear_transform(a, b, c, d);
    trans_sk_2 = s_k.fractional_linear_transform(a, b, c, d);
    trans_s_2 = s.fractional_linear_transform(a, b, c, d);
    [trans_ej_1, trans_ej_2] = trans_sj_2.get_endpoints();

    trans_points_2 = zeros(0, 1);

    new_points = point_line_bisector_intersection(trans_s_2, trans_ej_1, trans_sk_2, trans_sj_2, 1);
    trans_points_2 = [trans_points_2; new_points];

    new_points = point_line_bisector_intersection(trans_s_2, trans_ej_2, trans_sk_2, trans_sj_2, 3);
    trans_points_2 = [trans_points_2; new_points];

    points_2 = (d*trans_points_2 - b)./(-c*trans_points_2 + a);

    % Add everything together
    % TO-DO: Fix concat bug here
    points = unique([points; points_1; points_2]);
end

function points = point_point_bisector_intersections(s, p1, p2, g1, g2, ...
    R_for_g1, R_for_g2)
    % Construct the bisector of p1 and p2, which is the perpendicular
    % bisector of p1 and p2
    p1_to_p2 = GeodesicSegment(p1, p2);
    [a,b, c, d] = p1_to_p2.find_flt_to_imag_axis();
    [e1, e2] = p1_to_p2.fractional_linear_transform(a, b, c, d).get_endpoints();
    trans_s = s.fractional_linear_transform(a, b, c, d);

    radius = sqrt(imag(e1) * imag(e2));
    bisector = GeodesicSegment(-radius, radius);

    % Calculate intersection with s and bisector extended bi-infinitely
    points = bisector.intersections_with_geodesic(trans_s, true, false);

    points = (d * points - b)./(-c * points + a);

    if ~isempty(points)
        p = points(1);
        if ~(g1.get_region_of_point(p) == R_for_g1 && ...
             g2.get_region_of_point(p) == R_for_g2)
            points = [];
        end
    end      
end

function points = line_line_bisectors_intersection(s, g1, g2)
    points = zeros(0, 1);

    g2_center = g2.get_center_on_real_line();
    g2_radius = g2.get_radius_from_center();
    bisector_center_1 = g2_center + g2_radius;
    bisector_radius_1 = sqrt(bisector_center_1^2 + g2_radius^2 - g2_center^2);
    if imag(bisector_radius_1) == 0
            bisector_1 = GeodesicSegment(bisector_center_1 - bisector_radius_1, ...
                                         bisector_center_1 + bisector_radius_1);
            points_1  = bisector_1.intersections_with_geodesic(s, true, false);
            if ~isempty(points_1)
                if g1.get_region_of_point(points_1(1)) == 2 && ...
                   g2.get_region_of_point(points_1(1)) == 2
                    points = [points; points_1];
                end
            end
    end

    bisector_center_2 = g2_center - g2_radius;
    bisector_radius_2 = sqrt(bisector_center_2^2 + g2_radius^2 - g2_center^2);
    if imag(bisector_radius_2) == 0
        bisector_2 = GeodesicSegment(bisector_center_2 - bisector_radius_2, ...
                                     bisector_center_2 + bisector_radius_2);
        points_2 = bisector_2.intersections_with_geodesic(s, true, false);

        if ~isempty(points_2)
            if g1.get_region_of_point(points_2(1)) == 2 && ...
               g2.get_region_of_point(points_2(1)) == 2
                points = [points; points_2];
            end
        end
    end
end

function points = point_line_bisector_intersection(s, p, g, g2, R_for_g2)
    [e1, e2] = s.get_endpoints();
    min_x_on_s = min(real(e1), real(e2));
    max_x_on_s = max(real(e1), real(e2));

    u = real(p);
    v = imag(p);
    r = s.get_radius_from_center();
    c = s.get_center_on_real_line();

    C_1 = c - u;
    C_2 = r^2 - c^2 + u^2 - v^2;
    C_3 = u*v^2;
    
    A = C_1^2;
    B = C_1*C_2 - 2*C_3;
    C = C_2^2/4 + u*C_3;

    discriminant = B^2 - 4*A*C;
    points = [];
    if discriminant >= 0
        x_1 = (-B + sqrt(discriminant)) / (2*A);
        x_2 = (-B - sqrt(discriminant)) / (2*A);
        y_1_squared = r^2 - (x_1 - c)^2;
        y_2_squared = r^2 - (x_2 - c)^2;
        if y_1_squared > 0 && min_x_on_s <= x_1 && x_1 <= max_x_on_s
            candidate_point = x_1 + sqrt(y_1_squared)*1i;
            if g.get_region_of_point(candidate_point) == 2 && ...
               g2.get_region_of_point(candidate_point) == R_for_g2
                points = [points; candidate_point];
            end
        end
        if y_2_squared > 0 && min_x_on_s <= x_2 && x_2 <= max_x_on_s
            candidate_point = x_2 + sqrt(y_2_squared)*1i;
            if g.get_region_of_point(candidate_point) == 2 && ...
               g2.get_region_of_point(candidate_point) == R_for_g2
                points = [points; candidate_point];
            end
        end
    end
end