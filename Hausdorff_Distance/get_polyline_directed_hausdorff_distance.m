% calculate the directed Hausdorff distance
function polyline_directed_hausdorff_distance_polyline = get_polyline_directed_hausdorff_distance(true_geodesic,quasi_geodesic)
    % pruning
    polyline = get_unprunable_seg(true_geodesic,quasi_geodesic);
        
    % create segments on true geodesic from the intersections
    intersection_segments = get_voronoi_segments(true_geodesic, quasi_geodesic);

    % calculate the directed hausdorff distance between the segments on the true geodesic and their corresponding segment on pruned quasi-geodesic
    directed_hausdorff_distance = -1; % initialize
    for i = 1:length(intersections_segments)
        dist = get_geodesic_geodesic_directed_hausdorff(intersection_segments(i),polyline(i));
        if dist > directed_hausdorff_distance
            directed_hausdorff_distance = dist;
        end
    end
end
