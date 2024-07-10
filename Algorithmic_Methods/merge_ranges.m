% merge multiple ranges
function mergedRange = merge_ranges(ranges)
    % Initialize the intersection range
    mergedRange = [-inf, inf];
    
    % Find the intersection of all ranges
    for i = 1:size(ranges, 1)
        curr_range = ranges(i, :);
        if curr_range(1) > 2*pi
            curr_range = curr_range - [2*pi, 2*pi];
        end
        if curr_range(2) > 2*pi
            curr_range = [curr_range(1), 2*pi;
                          0, curr_range(2) - 2*pi];
        end
        mergedRange_size = size(mergedRange, 1);
        curr_range_size = size(curr_range, 1);
        new_mergedRange = zeros(0, 2);
        for j = 1:mergedRange_size
            for k = 1:curr_range_size
                r_lower = max(mergedRange(j, 1), curr_range(k, 1));
                r_upper = min(mergedRange(j, 2), curr_range(k, 2));
                if r_upper > r_lower
                    new_mergedRange = [new_mergedRange ; [r_lower, r_upper]];
                end
            end
        end
        mergedRange = new_mergedRange;
    end
    
    % Check if the intersection is valid
    %if isempty(mergedRange)
        %double(ranges)
        %error('The ranges do not overlap.');
    %end
end