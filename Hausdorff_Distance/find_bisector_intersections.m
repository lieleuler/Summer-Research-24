
% This function finds the intersections of a segment s with the bisector
% of two segments s_j and s_k. It first calculates all 14 possible
% intersections between s and the bisector, and then rules out candidates
% based on their coordinates
function points = find_bisector_intersections(s, s_j, s_k)

    % Step 1: Get candidate points from point-point bisectors. 
    candidate_points = zeros(0, 3);

    [ej_1, ej_2] = s_j.get_endpoints();
    [ek_1, ek_2] = s_k.get_endpoints();

    new_candidate_points = point_point_bisector_intersections(s, ej_1, ek_1);
    if ~isempty(candidate_points)
        candidate_points = [candidate_points; new_candidate_points(1), 1, 1];
    end

    new_candidate_points = point_point_bisector_intersections(s, ej_1, ek_2);
    if ~isempty(candidate_points)
        candidate_points = [candidate_points; new_candidate_points(1), 1, 3];
    end

    new_candidate_points = point_point_bisector_intersections(s, ej_2, ek_1);
    if ~isempty(candidate_points)
        candidate_points = [candidate_points; new_candidate_points(1), 3, 1];
    end

    new_candidate_points = point_point_bisector_intersections(s, ej_2, ek_2);
    if ~isempty(candidate_points)
        candidate_points = [candidate_points; new_candidate_points(1), 3, 3];
    end

    % Step 2: Get candidate points from line-line bisectors
    new_candidate_points = line_line_bisectors_intersection(s, s_j, s_k);
    for p = new_candidate_points % TO-DO: Make this not a for loop!!!
        candidate_points = [candidate_points; p, 2, 2];
    end

end

function points = point_point_bisector_intersections(s, p1, p2)
    % Construct the bisector of p1 and p2, which is the perpendicular
    % bisector of p1 and p2
    p1_to_p2 = GeodesicSegment(p1, p2);
    half_length = p1_to_p2.get_length() / 2;
    midpoint = p1_to_p2.travel_from_start(half_length);
    perpendicular_dir = pi / 2 + ...
                        p1_to_p2.get_angle_with_vertical(half_length);
    bisector = GeodesicSegment(midpoint, ...
                               get_point_along_direction(midpoint, ...
                                                         perpendicular_dir, ...
                                                         1));

    % Calculate intersection with s and bisector extended bi-infinitely
    points = bisector.intersections_with_geodesic(s, true, false);
end

function points = line_line_bisectors_intersection(s, g1, g2)
    [a, b, c, d] = g1.find_flt_to_imag_axis();
    trans_g2 = g2.fractional_linear_transform(a, b, c, d);
    trans_s = s.fractional_linear_transform(a, b, c, d);

    trans_g2_center = trans_g2.get_center_on_real_line();
    trans_g2_radius = trans_g2.get_radius_from_center();
    bisector_center_1 = trans_g2_center + trans_g2_radius;
    bisector_radius_1 = sqrt(bisector_center_1^2 + trans_g2_radius^2 - trans_g2_center^2);
    bisector_center_2 = trans_g2_center - trans_g2_radius;
    bisector_radius_2 = sqrt(bisector_center_2^2 + trans_g2_radius^2 - trans_g2_center^2);
    bisector_1 = GeodesicSegment(bisector_center_1 - bisector_radius_1, ...
                                 bisector_center_1 + bisector_radius_1);
    bisector_2 = GeodesicSegment(bisector_center_2 - bisector_radius_2, ...
                                 bisector_center_2 + bisector_radius_2);

    points_1  = bisector_1.intersections_with_geodesic(trans_s, true, false);
    points_2 = bisector_2.intersections_with_geodesic(trans_s, true, false);
    points = [points_1; points_2];
    points = (d * points - b)./(-c*points + a);
end

function points = line_point_bisector_intersection(s, p, g)
    [a, b, c, d] = g.find_flt_to_imag_axis();
    trans_p = (a*p + b) / (c*p + d);
    trans_s = s.fractional_linear_transform(a, b, c, d);

    % ===================== %
    % == INSERT EQUATION == %
    % ===================== %

    u = real(trans_p);
    v = imag(trans_p);
    r = trans_s.get_radius_from_center();
    c = trans_s.get_center_on_real_line();

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
        y_1_squared = sqrt(r^2 - (x_1 - c)^2);
        y_2_squared = sqrt(r^2 - (x_2 - c)^2);
        if y_1_squared > 0
            points = [points; x_1 + sqrt(y_1_squared)*1i];
        end
        if y_2_squared > 0
            points = [points; x_2 + sqrt(y_2_squared)*1i];
        end
    end
end