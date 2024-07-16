% constants for quasi-geodesic
% epsilon varies
DELTA = log(1 + sqrt(2));
lambda = 10;
target_length = 3; 
step_size = 0.1;
min_segment_splits = 1;
pruning_sample_ratio = 100; % 1/N of the segments will be sampled for pruning

iteration_num = 40;
dist_hausdorff = [];
length_quasi_geodesic = [];

% variable
epsilon_seq = 2:4:34;

% data storage
hausdorff_distance_data = [];
hausdorff_distance_data_mean = [];
hausdorff_distance_data_median = [];
hausdorff_distance_data_min = [];
hausdorff_distance_data_max = [];

for eps = epsilon_seq
    lap_hausdorff_distance_data = [];

    for i = 1:iteration_num
        disp("Quasi-geodesic NO."+i+" @ epsilon = "+eps);
    
        tic
        [start_points, end_points] = random_walk_hyperbolic(lambda, eps, step_size, target_length, min_segment_splits,0); % generate quasi-geodesic
        quasi_geodesic = points_to_quasigeodesic(start_points, end_points);
        num_points = length(quasi_geodesic);
        pruning_sample = ceil(num_points/pruning_sample_ratio);
        true_geodesic = GeodesicSegment(start_points(1),end_points(num_points));
        verification = verify(quasi_geodesic, true_geodesic, DELTA, lambda, eps, pruning_sample); % verify the hausdorff distance
        toc
    
        % store data
        lap_hausdorff_distance_data = [lap_hausdorff_distance_data, verification(2)];
        
        disp("Length of geodesic = "+target_length);
        disp("Number of steps taken = "+num_points);
        disp("R bound = "+verification(1)+"; Hausdorff Distance = "+verification(2));
        disp("--------------------------");
        disp(' ');
        disp(' ');
        
    end
    hausdorff_distance_data = [hausdorff_distance_data,lap_hausdorff_distance_data];
    hausdorff_distance_data_mean = [hausdorff_distance_data_mean, mean(lap_hausdorff_distance_data)];
    hausdorff_distance_data_median = [hausdorff_distance_data_median, median(lap_hausdorff_distance_data)];
    hausdorff_distance_data_min = [hausdorff_distance_data_min, min(lap_hausdorff_distance_data)];
    hausdorff_distance_data_max = [hausdorff_distance_data_max, max(lap_hausdorff_distance_data)];
end



% Store the values in a structure
results.hausdorff_distance_data = hausdorff_distance_data;
results.hausdorff_distance_data_mean = hausdorff_distance_data_mean;
results.hausdorff_distance_data_median = hausdorff_distance_data_median;
results.hausdorff_distance_data_min = hausdorff_distance_data_min;
results.hausdorff_distance_data_max = hausdorff_distance_data_max;
results.epsilon_seq = epsilon_seq;
results.iteration_size = iteration_num;


% Save the results to a .mat file
folder_path = '/Users/leesmackintosh/Desktop/Summer Research/Summer-Research-24-main/Results';
save(fullfile(folder_path, 'epsilon_test_results.mat'), '-struct', 'results');
