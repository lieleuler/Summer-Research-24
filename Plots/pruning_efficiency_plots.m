% Load the results from the .mat file in the specified folder
folder_path = '/Users/leesmackintosh/Desktop/Summer Research/Summer-Research-24-main/Results'; % Change this to your desired folder path
loaded_results = load(fullfile(folder_path, 'pruning_efficiency_test_results.mat'));

% plot computational requirement
figure;

lower_error = loaded_results.computation_requirement_mean - loaded_results.computation_requirement_min;
upper_error = loaded_results.computation_requirement_max - loaded_results.computation_requirement_mean;

errorbar(loaded_results.pruning, loaded_results.computation_requirement_mean, lower_error, upper_error, 's', ...
    'MarkerSize', 10, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red','Color','black','LineWidth', 2);


xlabel('Pruning Ratio = 1/N');
ylabel('Computation Requirement');
title('Pruning Efficiency: Computation Requirement');
grid on;
hold off;

% plot run time
figure;

lower_error = loaded_results.runTime_mean - loaded_results.runTime_min;
upper_error = loaded_results.runTime_max - loaded_results.runTime_mean;

errorbar(loaded_results.pruning, loaded_results.runTime_mean, lower_error, upper_error, 's', ...
    'MarkerSize', 10, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red','Color','black','LineWidth', 2);


xlabel('Pruning Ratio = 1/N');
ylabel('Run Time');
title('Pruning Efficiency: Run Time');
grid on;
hold off;