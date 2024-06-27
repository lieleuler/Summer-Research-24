
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
    transformed_g2 = g2.fractional_linear_transform(a, b, c, d);
    transformed_s = s.fractional_linear_transform(a, b, c, d);

    transformed_g2_center = transformed_g2.get_center_on_real_line();
    transformed_g2_radius = transformed_g2.get_radius_from_center();
    bisector_center_1 = transformed_g2_center + transformed_g2_radius;
    bisector_radius_1 = sqrt(bisector_center_1^2 + transformed_g2_radius^2 - transformed_g2_center^2);
    bisector_center_2 = transformed_g2_center - transformed_g2_radius;
    bisector_radius_2 = sqrt(bisector_center_2^2 + transformed_g2_radius^2 - transformed_g2_center^2);
    bisector_1 = GeodesicSegment(bisector_center_1 - bisector_radius_1, ...
                                 bisector_center_1 + bisector_radius_1);
    bisector_2 = GeodesicSegment(bisector_center_2 - bisector_radius_2, ...
                                 bisector_center_2 + bisector_radius_2);

    points_1  = bisector_1.intersections_with_geodesic(transformed_s, true, false);
    points_2 = bisector_2.intersections_with_geodesic(transformed_s, true, false);
    points = [points_1; points_2];
    points = (d * points - b)./(-c*points + a);
end

function line_point_bisector_intersection(s, p, g)
    [a, b, c, d] = g.find_flt_to_imag_axis();
    trans_p = (a*p + b) / (c*p + d);
    trans_s = s.fractional_linear_transform(a, b, c, d);

    % ===================== %
    % == INSERT EQUATION == %
    % ===================== %
end