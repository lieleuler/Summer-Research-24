
function point = generateJumpPoint(z, jump_ranges_in_eps_circle, total_weight)
    % Pick a random portion of the range by weight
    i = 0;
    random_value = total_weight * rand();
    while true
        i = i + 1;

        curr_weight = jump_ranges_in_eps_circle(i, 4);
        random_value = random_value - curr_weight;
        if random_value <= 0
            break
        end
    end
    % Pick random point on that range
    theta_1 = jump_ranges_in_eps_circle(i, 1);
    theta_2 = jump_ranges_in_eps_circle(i, 2);
    radius = jump_ranges_in_eps_circle(i, 3);

    theta = theta_1 + rand()*(theta_2 - theta_1);

    point = get_point_along_direction(z, theta, radius);
end