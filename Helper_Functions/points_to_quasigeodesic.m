

function quasi_geodesic = points_to_geodesic(points)
    n = length(points)-1;
    quasi_geodesic(1, n) = GeodesicSegment;
    for i = 1:n
        quasi_geodesic(i) = GeodesicSegment(points(i), points(i + 1));
    end
end
