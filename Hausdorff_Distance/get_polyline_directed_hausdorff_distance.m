% calculate the directed Hausdorff distance
function [polyline_directed_hausdorff_distance, intersections_segments, corresponding_segments] = get_polyline_directed_hausdorff_distance(true_geodesic, quasi_geodesic)
    % pruning
    polyline = get_unprunable_seg(true_geodesic,quasi_geodesic, 5);
        
    % create segments on true geodesic from the intersections
    [intersections_segments, corresponding_segments] = get_voronoi_segments(true_geodesic, quasi_geodesic);

    % calculate the directed hausdorff distance between the segments on the true geodesic and their corresponding segment on pruned quasi-geodesic
    polyline_directed_hausdorff_distance = -1; % initialize

    best_segment = GeodesicSegment();
    best_intersection = 0;
    for i = 1:length(intersections_segments)
        dist = get_geodesic_directed_hausdorff(intersections_segments(i), corresponding_segments(i));
        if dist > polyline_directed_hausdorff_distance
            best_segment = corresponding_segments(i);
            best_intersection = intersections_segments(i);
            polyline_directed_hausdorff_distance = dist;
        end
    end
end
