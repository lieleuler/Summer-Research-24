% constants for quasi-geodesic
% pruning_sample_ratio varies
DELTA = log(1 + sqrt(2));
eps = 5;
lambda = 10;
target_length = 2; 
step_size = 0.1;
min_segment_splits = 1;

iteration_num = 5;
dist_hausdorff = [];
length_quasi_geodesic = [];

% variable
pruning = 2:4:34;

% data storage
runTime = [];
computation_requirement = [];
computation_requirement_mean = [];
computation_requirement_median = [];
computation_requirement_min = [];
computation_requirement_max = [];

runTime_mean = [];
runTime_median = [];
runTime_min = [];
runTime_max = [];


for pruning_sample_ratio = pruning % 1/N of the segments will be sampled for pruning
    
    lap_runTime = [];
    lap_computation_requirement = [];
    
    for i = 1:iteration_num
        disp("Quasi-geodesic NO."+i+" @ pruning sample ratio = "+pruning_sample_ratio);
    
        startTime = cputime;
        [start_points, end_points] = random_walk_hyperbolic(lambda, eps, step_size, target_length, min_segment_splits,0); % generate quasi-geodesic
        quasi_geodesic = points_to_quasigeodesic(start_points, end_points);
        num_points = length(quasi_geodesic);
        pruning_sample = ceil(num_points/pruning_sample_ratio);
        true_geodesic = GeodesicSegment(start_points(1),end_points(num_points));
        verification = verify(quasi_geodesic, true_geodesic, DELTA, lambda, eps, pruning_sample); % verify the hausdorff distance
        elapsedTime = cputime - startTime;
    
        % store data
        lap_runTime = [lap_runTime, elapsedTime];
        lap_computation_requirement = [lap_computation_requirement, num_points * pruning_sample];
        
        disp("Length of geodesic = "+target_length);
        disp("Number of steps taken = "+num_points);
        disp("R bound = "+verification(1)+"; Hausdorff Distance = "+verification(2));
        disp("Elapsed Time: "+elapsedTime);
        disp("Computation Requirement: "+num_points * pruning_sample);
        disp("--------------------------");
        disp(' ');
        disp(' '); 
    end

    % stat
    computation_requirement_mean = [computation_requirement_mean, mean(lap_computation_requirement)];
    computation_requirement_median = [computation_requirement_median, median(lap_computation_requirement)];
    computation_requirement_min = [computation_requirement_min, min(lap_computation_requirement)];
    computation_requirement_max = [computation_requirement_max, max(lap_computation_requirement)];

    runTime_mean = [runTime_mean, mean(lap_runTime)];
    runTime_median = [runTime_median, median(lap_runTime)];
    runTime_min = [runTime_min, min(lap_runTime)];
    runTime_max = [runTime_max, max(lap_runTime)];
end

% Store the values in a structure
results.pruning = pruning;
results.computation_requirement_mean = computation_requirement_mean;
results.computation_requirement_median = computation_requirement_median;
results.computation_requirement_min = computation_requirement_min;
results.computation_requirement_max = computation_requirement_max;

results.runTime_mean = runTime_mean;
results.runTime_median = runTime_median;
results.runTime_min = runTime_min;
results.runTime_max = runTime_max;


% Save the results to a .mat file
folder_path = '/Users/leesmackintosh/Desktop/Summer Research/Summer-Research-24-main/Results';
save(fullfile(folder_path,'pruning_efficiency_test_results.mat'), '-struct', 'results');



