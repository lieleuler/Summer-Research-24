% Load the results from the .mat file in the specified folder
folder_path = '/Users/leesmackintosh/Desktop/Summer Research/Summer-Research-24-main/Results'; % Change this to your desired folder path
loaded_results = load(fullfile(folder_path, 'lambda_test_results.mat'));

disp("Iteration Size: "+ loaded_results.iteration_size);

% plot
figure;

lower_error = loaded_results.hausdorff_distance_data_mean - loaded_results.hausdorff_distance_data_min;
upper_error = loaded_results.hausdorff_distance_data_max - loaded_results.hausdorff_distance_data_mean;

hold on

hMean = errorbar(loaded_results.lambda_seq, loaded_results.hausdorff_distance_data_mean, lower_error, upper_error, 's', ...
    'MarkerSize', 10, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red','Color','black','LineWidth', 2);

hMedian = plot(loaded_results.lambda_seq, loaded_results.hausdorff_distance_data_median, 'p',...
    'MarkerSize', 12, 'MarkerEdgeColor', 'green', 'MarkerFaceColor', 'green');

xlabel('Lambda');
ylabel('Hausdorff Distance');
title('Hausdorff Distance with Various Lambda');
grid on;
legend([hMean, hMedian], {'Mean', 'Median'}, 'Location', 'best');
hold off;