
% calculate the intersection of voronoi diagram with the true geodesic, and
% split intersections into segments
function split_segment = get_voronoi_segments(geodesic, pruned_polyline)
    first_seg = pruned_polyline(1);
    split_segment_map = [geodesic, first_seg];
    for i = 1:size(pruned_polyline, 1)
        new_polyline_seg = pruned_polyline(i);
        new_split_segment_map = zeros(0, 1);
        for j = 1:size(split_segment_map, 1)
            old_polyline_seg = split_segment_map(j, 2);
            corresponding_geodesic_seg = split_segment_map(j, 1);
            intersections = find_bisector_intersections(corresponding_geodesic_seg, ...
                                                        old_polyline_seg, ...
                                                        new_polyline_seg);
            new_split_segment_map = [new_split_segment_map; 
                                 split_segment_with_intersections(corresponding_geodesic_seg, ...
                                                                  intersections, ...
                                                                  new_polyline_seg, ...
                                                                  old_polyline_seg)];
        end
        split_segment_map = new_split_segment_map;
    end
    split_segment = split_segment_map(:, 1);
end

function split_segment_map = split_segment_with_intersections(segment, intersections, new_s, old_s)
    if isempty(intersections)
        if new_s.dist_from_point(sample_point) <= old_s.dist_from_point(sample_point)
            split_segment_map = [segment, new_s];
        else
            split_segment_map = [segment, old_s];
        end
        return
    end

    % Sort intersections in order of 
    intersections = sort(intersections, "ComparisonMethod", "real");

    end_points = segment.get_endpoints(); % get end points of the true geodesic
    start_point = end_points(1); % initialize

    intersection_count = length(intersections);
    split_segment_map = zeros(intersection_count, 2);
    for i = 1:(intersection_count + 1)
        if i <= intersection_count
            end_point = intersections(i);
        else
            end_point = end_points(2);
        end
        % connecting each intersections on the true geodesic
        new_subsegment = GeodesicSegment(start_point, end_point);
        split_segment_map(i, 1) = new_subsegment;
        if new_s.dist_from_point(sample_point) <= old_s.dist_from_point(sample_point)
            split_segment_map(i, 2) = new_s;
        else
            split_segment_map(i, 2) = old_s;
        end

        % Update start point to be the previous endpoint
        start_point = end_point;
    end
end
