
function quasi_geodesic = points_to_quasigeodesic(start_points, end_points)
    n = length(start_points);
    quasi_geodesic(1, n) = GeodesicSegment;
    for i = 1:n
        quasi_geodesic(i) = GeodesicSegment(start_points(i), end_points(i));
    end
end
