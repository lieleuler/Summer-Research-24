% generate the Voronoi diagram for quasi_geodesic polyline
function voronoi = get_voronoi(quasi_geodesic)
    %TODO

end

% calculate the intersection of voronoi diagram with the true geodesic
function voronoi_intersection = get_voronoi_intersection(voronoi,geodesic)
    %TODO

end


% calculate the minimal distance MD(geo_1,geo_2) for pruning
function minimal_distance = get_minimal_distance(geodesic_1,geodesic_2,c,r) 
    % mobius transformation for geodesic_1
    geodesic_1 = mobius_transform(geodesic_1); % transform geodesic_1 to a straight line on imaginary axis
    
    c = center(geodesic_2); % calculate the center of the circle on which geodesic_2 lies
    r = radius(geodesic_2); % calculate the radius of the circle on which geodesic_2 lies

    syms si % angle at center of d2 from x-axis to the chosen point p on s

    % calculate the range of si based on end points of the segement
    start_pt_1 = get_start_point(geodesic_1);
    end_pt_1 = get_end_point(geodesic_1);
    
    range_si = sort(atan(-start_pt_(2)/c),atan(-end_pt_1(2)/c));

    % calculate the range of theta based on end points of the segement
    start_pt_2 = get_start_point(geodesic_2);
    end_pt_2 = get_end_point(geodesic_2);

    range_theta = sort(atan((start_pt_2(1)-c)/start_pt_2(2)),atan((end_pt_2(1)-c)/end_pt_2(2)));

    % take derivative of distance function
    % Expression for x = cos(theta)
    x = (2*r*c*cos(si)^2) / (c^2*cos(si)^2 + r^2);
    % Derivative of x with respect to si
    dx_dsi = diff(x, si);

    % Substitute x and x' in the given expression
    combined_expr = (1/(2*c*r*(1 - x^2)^(3/2))) * csc(2*si)^2 * (4*c*(c - 2*r*x)*(1 - x^2) ...
        + 4*(c^2 + 2*r^2 - 2*c*r*x)*(1 - x^2)*cos(2*si) ...
        + 2*(-2*c*r + c^2*x + 2*r^2*x)*dx_dsi*sin(2*si) ...
        + c*(-2*r + c*x)*dx_dsi*sin(4*si));

    % Simplify the expression
    combined_expr = simplify(combined_expr);

    % Define the equation to solve
    eqn = combined_expr == 0;

    % Use numerical solver to find si in the range range_si
    si_roots = vpasolve(eqn, si, range_si);
    x_si = (2*r*c*cos(si_roots)^2) / (c^2*cos(si_roots)^2 + r^2); % evaluate numerical value of x = cos(theta) at si_roots
    theta_si = acos(x_si); % evaluate numerical value of theta at si_roots
    
    % output MD
    if isempty(si_roots) % no critical point in the specified range_si
        minimal_distance = min(dist_H(start_pt_1,start_pt_2),dist_H(start_pt_1,end_pt_2),dist_H(start_pt_2,end_pt_1),dist_H(end_pt_1,end_pt_2));
    elseif theta_si > range_theta(2) || theta_si < range_theta(1) % theta_si outside the range_theta 
        minimal_distance = min(dist_H(start_pt_1,start_pt_2),dist_H(start_pt_1,end_pt_2),dist_H(start_pt_2,end_pt_1),dist_H(end_pt_1,end_pt_2));
    else % normal case 
        minimal_distance = acosh(-1+(2*(c^2 - 2*c*r*x_si) * (cos(si_roots))^2 + 2*r^2)/(c*r * sqrt(1-x_si^2) * sin(2*si_roots)));
    end
end

% =========================================================================================================

% sample input
sampling_size = 100; % number of sample to find the closeSeg
quasi_geodesic = [geodesic_1, geodesic_2, geodesic_3]; % quasi-geodesic consists of geodesic segments
true_geodesic = geodesic; % the true geodesic connecting the end points of the quasi-geodesic
num_step = length(quasi_geodesic); % number of segments in our quasi-geodesic equals the number of steps we ran


% pruning method for the directed Hausdorff distance algorithm
function unprunable_seg = get_unprunable_seg(true_geodesic,quasi_geodesic)
    % uniformly sample geodesic segments and find minimal dist2closeSeg
    dist2closeSeg = infinity; % initialization
    closeSeg_index = -1; % initialization
    
    assert(sampling_size < num_step, "sampling size too big.");

    for i = 1:sampling_size:num_step
        dist_hausdorff = get_geodesic2geodesic_directed_hausdorff(true_geodesic,quasi_geodesic(i));
        if dist_hausdorff < dist2closeSeg
            dist2closeSeg = dist_hausdorff;
            closeSeg_index = i;
        end
    end

    assert(closeSeg_index > 0, "can't find a closeSeg. something is wrong.");

    % iterate through all geodesic segments to compare their MD with dist2closeSeg for pruning
    unprunable_seg = [];
    unprunable_seg_index = []; % index for segments not pruned
    for i = 1:num_step
        if get_minimal_distance(true_geodesic, quasi_geodesic(i)) <= dist2closeSeg % not safe for pruning
            unprunable_seg_index = [unprunable_seg_index,i];
            unprunable_seg = [unprunable_seg,quasi_geodesic(i)];
        end
    end
    
    assert(isempty(unprunable_seg), "everything is pruned.");
end

% calculate the directed Hausdorff distance
function directed_hausdorff_distance = get_directed_hausdorff_distance(true_geodesic,quasi_geodesic)
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
    directed_hausdorff_distance = -1; % initialize
    for i = 1:length(intersections_seg_index)
        dist = get_geodesic_geodesic_directed_hausdorff(intersection_segments(i),polyline(i));
        if dist > directed_hausdorff_distance
            directed_hausdorff_distance = dist;
        end
    end

end
