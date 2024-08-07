
% calculate the intersection of voronoi diagram with the true geodesic, and
% split intersections into segments
function [split_segments, corresponding_segments] = get_voronoi_segments(geodesic, pruned_polyline)
    first_seg = pruned_polyline(1);
    split_segment_map = [geodesic, first_seg];
    for i = 2:length(pruned_polyline)
        new_polyline_seg = pruned_polyline(i);
        new_split_segment_map = zeros(0, 2);
        for j = 1:height(split_segment_map)
            corresponding_geodesic_seg = split_segment_map(j, 1);
            old_polyline_seg = split_segment_map(j, 2);
            intersections = find_bisector_intersections( ...
                                corresponding_geodesic_seg, ...
                                old_polyline_seg, ...
                                new_polyline_seg);
            new_split_segment_map = [new_split_segment_map; 
                                     split_segment_with_intersections( ...
                                                corresponding_geodesic_seg, ...
                                                intersections, ...
                                                new_polyline_seg, ...
                                                old_polyline_seg)];
        end
        split_segment_map = new_split_segment_map;

        % === TESTING === %
        %colors = ["r", "g", "b", "c", "m", "y"];
        %colors_len = length(colors);
        
        %hold on
        %for k = 1:height(split_segment_map)
        %    color = colors(mod(k - 1, colors_len) + 1);
        %    split_segment_map(k, 1).plot(100, color)
        %    split_segment_map(k, 2).plot(100, color)
        %    GeodesicSegment(split_segment_map(k, 1).get_midpoint(), ...
        %                    split_segment_map(k, 2).get_midpoint()).plot(100, "k")
        %end
        %"H"
        %waitforbuttonpress
        %clf
        % =============== %
        % Bugs:
        %  - First 3 [X]
        %  - Yellow Segment [X]
        %  - Point/Line Voronoi? [X]

    end
    split_segments = split_segment_map(:, 1);
    corresponding_segments = split_segment_map(:, 2);
end

function split_segment_map = split_segment_with_intersections(segment, intersections, new_s, old_s)
    % Case with no intersections
    if isempty(intersections)
        sample_point = segment.get_midpoint();
        if old_s.dist_from_point(sample_point) <= new_s.dist_from_point(sample_point)
            split_segment_map = [segment, old_s];
        else
            split_segment_map = [segment, new_s];
        end
        return
    end

    % TO-DO: Write this comment
    intersections = sort(unique(intersections), "ComparisonMethod", "real");

    [seg_start_point, seg_end_point] = segment.get_endpoints(); % get end points of the true geodesic
    start_point = seg_start_point; % initialize

    intersection_count = length(intersections);
    split_segment_map(intersection_count + 1, 2) = GeodesicSegment;
    for i = 1:(intersection_count + 1)
        if i <= intersection_count
            end_point = intersections(i);
        else
            end_point = seg_end_point;
        end

        % Connecting each intersection on the true geodesic
        new_subsegment = GeodesicSegment(start_point, end_point);
        split_segment_map(i, 1) = new_subsegment;
        sample_point = new_subsegment.get_midpoint();
        if old_s.dist_from_point(sample_point) <= new_s.dist_from_point(sample_point)
            split_segment_map(i, 2) = old_s;
        else
            split_segment_map(i, 2) = new_s;
        end

        % Update start point to be the previous endpoint
        start_point = end_point;
    end
end
