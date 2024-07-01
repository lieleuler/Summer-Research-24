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
