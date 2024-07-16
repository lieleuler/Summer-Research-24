% Load the results from the .mat file in the specified folder
folder_path = '/Users/leesmackintosh/Desktop/Summer Research/Summer-Research-24-main/Results'; % Change this to your desired folder path
loaded_results = load(fullfile(folder_path, 'target_length_test_results.mat'));

disp("Iteration Size: "+ loaded_results.iteration_size);

% plot
figure;

lower_error = loaded_results.hausdorff_distance_data_mean - hausdorff_distance_data_min;
upper_error = loaded_results.hausdorff_distance_data_max - hausdorff_distance_data_mean;

hold on

hMean = errorbar(loaded_results.target_length_seq, loaded_results.hausdorff_distance_data_mean, lower_error, upper_error, 's', ...
    'MarkerSize', 10, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red','Color','black','LineWidth', 2);

hMedian = plot(loaded_results.target_length_seq, loaded_results.hausdorff_distance_data_mean, 'p',...
    'MarkerSize', 12, 'MarkerEdgeColor', 'green', 'MarkerFaceColor', 'green');

xlabel('Target Length');
ylabel('Hausdorff Distance');
title('Hausdorff Distance with Various Target Length');
grid on;
legend([hMean, hMedian], {'Mean', 'Median'}, 'Location', 'best');
hold off;