% generate branching_num for random tree walk algorithm
function branching_num = branching_z_function(starting_branch_num, steepness, midpoint, step_size, goal_branches)
    A = starting_branch_num; % starting y-value
    B = 1; % ending y-value
    C = midpoint; % middle point of the slope
    k = steepness; % slope
    branching_num = [];
    total_nodes = 0;
    for x = 1:step_size
        f(x) = ceil(A + (B-A)/(1+exp(-k*(x-C))));
        branching_num = [branching_num, f(x)];
        total_nodes = total_nodes + f(x);
    end
end
