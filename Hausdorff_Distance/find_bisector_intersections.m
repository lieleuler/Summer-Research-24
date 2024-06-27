
% This function finds the intersections of a segment s with the bisector
% of two segments s_j and s_k. It first calculates all 14 possible
% intersections between s and the bisector, and then rules out candidates
% based on their coordinates
function points = find_bisector_intersections(s, s_j, s_k)
    % Step 1: Get candidate points from point-point bisectors. 
    points = zeros(0, 1);

    [ej_1, ej_2] = s_j.get_endpoints();
    [ek_1, ek_2] = s_k.get_endpoints();

    candidate_points = point_point_bisector_intersections(s, ej_1, ek_1);
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