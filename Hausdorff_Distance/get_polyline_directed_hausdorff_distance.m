% calculate the directed Hausdorff distance
function polyline_directed_hausdorff_distance_polyline = get_polyline_directed_hausdorff_distance(true_geodesic,quasi_geodesic)
    % pruning
    polyline = get_unprunable_seg(true_geodesic,quasi_geodesic);
    
    % draw the voronoi diagram
    voronoi = get_voronoi(polyline);

    % calculate the intersection of voronoi diagram with the true geodesic
    intersections_info = get_voronoi_intersection(voronoi,true_geodesic); % assumed output: [[inter_1,seg_i],[inter_2,seg_j],...]
    intersections = cellfun(@(c) c{1}, intersections_info); % extract intersection only
    intersections_seg_index = cellfun(@(c) c{2}, intersections_info); % extract index of corresponding segments of voronoi cell only
    
    % create segments on true geodesic from the intersections
    intersection_segments = [];
    end_points = get_end_points_from_geodesic(true_geodesic); % get end points of the true geodesic
    start_point = end_points(1); % initialize
    for i = 1:length(intersections)
        % connecting each intersections on the true geodesic
        intersection_segments = [intersection_segments, get_geodesic_from_end_points(start_point,intersections(i))];
        start_point = intersections(i);
    end
    end_point = end_points(2);
    intersection_segments = [intersection_segments, get_geodesic_from_end_points(intersections(end),end_point)];

    % calculate the directed hausdorff distance between the segments on the true geodesic and their corresponding segment on pruned quasi-geodesic
    polyline_directed_hausdorff_distance = -1; % initialize
    for i = 1:length(intersections_seg_index)
        dist = get_geodesic_directed_hausdorff(intersection_segments(i),polyline(i));
        if dist > polyline_directed_hausdorff_distance
            polyline_directed_hausdorff_distance = dist;
        end
    end
end
