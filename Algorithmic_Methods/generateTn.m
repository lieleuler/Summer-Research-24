% generate bounded t_n using randomization
function [tn, phi] = generateTn(t,range,step_size) % ?might be able to generate "geodesic"
    % Generate random value within size of range for uniform distribution
    range_size = 0;
    for i = 1:size(range, 1)
        range_size = range_size + (range(i, 2) - range(i, 1));
    end

    % Get phi via random value
    random_value = range_size * rand();
    for i = 1:size(range, 1)
        random_value = random_value - (range(i, 2) - range(i, 1));
        if random_value <= 0
            phi = range(i, 2) + random_value;
            break
        end
    end

    % Calculate the coordinates of the random point
    tn = get_point_along_direction(t, phi, step_size);
end