function random_walk_hyperbolic(n_steps, step_size)
    % This function performs a random walk in the hyperbolic upper half-plane
    % n_steps: Number of steps in the random walk
    % step_size: Fixed hyperbolic distance to move at each step

    % Initialize the starting point
    z = 1i;  % Start at (0, 1) to avoid the real axis
    
    % Pre-allocate array for plotting
    points = zeros(n_steps + 1, 1);
    points(1) = z;
    
    % Perform the random walk
    for k = 1:n_steps
        % Get current point
        z = points(k);

        % Random angle in [0, 2*pi]
        theta = 2 * pi * rand();
        
        % Destination point on the geodesic
        % Moving along a circle arc centered at real axis
        t = linspace(0, step_size, 100);  % parameterize the distance along geodesic
        % Compute geodesic from z along angle theta
        u = real(z);
        v = imag(z);
        center_on_real_axis = v*(pi/2 - tan(theta)) + u; 
        geodesic = GeodesicSegment(z, z + 2*(center_on_real_axis - u));
        % Step along geodesic by the step_size
        if 0 <= theta && theta <= pi
            new_z = geodesic.get(step_size);
        else
            new_z = geodesic.get(-step_size);
        end
        points(k+1) = new_z;
    end
    
    % Plot the path
    figure;
    plot(real(points), imag(points), 'o-');
    xlabel('Real part');
    ylabel('Imaginary part');
    title('Random Walk in the Hyperbolic Upper Half Plane');
    axis equal;
end