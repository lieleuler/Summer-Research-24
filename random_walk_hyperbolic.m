
function random_walk_hyperbolic(n_steps, step_size)
    % This function performs a random walk in the hyperbolic upper half-plane
    % n_steps: Number of steps in the random walk
    % step_size: Fixed hyperbolic distance to move at each step

    lambda = 1.5;
    eps = 1;

    angle_divisions = 360;

    % Initialize the starting point
    z = 1i;  % Start at (0, 1) to avoid the real axis
    
    % Pre-allocate array for plotting
    points = zeros(n_steps + 1, 1);
    points(1) = z;
    
    % Perform the random walk
    for k = 1:n_steps
        % Get current point and the possible points to step to in a circle
        % around it
        z = points(k);
        points_on_step_circle = get_points_on_circle(z, step_size, angle_divisions);

        % Determine the directions we can go without violating
        % quasi-geodesic lower bound condition
        for t_i = 1:(k - 2)
            % Get the i-th previous segment and determine if 
            segment_i = GeodesicSegment(points(t_i), points(t_i + 1));

            % Determine the s-neighborhood around the segment, which is a
            % neighborhood within which we cannot guarentee stepping into
            % will preserve the quasi-geodesic lower bound
            s_neighborhood = get_S_ngbh(segment_i, t_i*step_size, ...
                k*step_size, lambda, eps, step_size);
            
            [point1, point2] = intersection_of_circle_and_s_ngbh( ...
                points_on_step_circle, s_neighborhood, angle_divisions);

            %% ADD CODE HERE FOR FINDING RANGE BETWEEN POINT 1 AND POINT 2 HERE %%
        end

        % Get random angle in [0, 2*pi]
        % TO-DO: Modify this code to get random angle in allowed range
        % (once the code for figuring out allowed range is written above)
        theta = 2 * pi * rand();
        
        % Compute geodesic from z in direction indicated theta
        points(k+1) = get_point_along_direction(z, theta, step_size);
    end
    
    % Plot the path
    figure;
    plot(real(points), imag(points), 'o-');
    xlabel('Real part');
    ylabel('Imaginary part');
    title('Random Walk in the Hyperbolic Upper Half Plane');
    axis equal;
end

function [int1, int2] = intersection_of_circle_and_s_ngbh(circle, ...
    s_neighborhood, angle_divisions)
    % Find a point inside the s-neighborhood. If one is found,
    % binary search to find the other boundary point
    t_after_boundary = -1;
    initial_point_in_ngbh = s_neighborhood(circle(1));
    for t = 2:angle_divisions
        if s_neighborhood(circle(t)) ~= initial_point_in_ngbh
            t_after_boundary = t;
            break
        end
    end
    % (Finding 2nd boundary point)
    if t_after_boundary ~= -1
        lower_bound = t_after_boundary;
        upper_bound = angle_divisions;
        while lower_bound < upper_bound
            mid = floor((lower_bound + upper_bound) / 2);
            if s_neighborhood(circle(t)) == initial_point_in_ngbh
                upper_bound = mid;
            else
                lower_bound = mid;
            end
        end
        if initial_point_in_ngbh
            if mid == angle_divisions
                mid = 0
            end
            int1 = circle(mid + 1);
            int2 = circle(t_after_boundary - 1);
        else 
            int1 = circle(t_after_boundary);
            int2 = circle(mid);
        end
    % This else case occurs when no boundary was found. In this case, the
    % circle does not intersect with the s-neighborhood, so we return [-1,
    % -1] to signify this
    else
        if initial_point_in_ngbh
            throw(MException("Circle entirely inside s-neighborhood!"))
        end
        int1 = -1;
        int2 = -1;
    end
end