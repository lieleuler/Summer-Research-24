
function random_walk_hyperbolic(n_steps, step_size)
    % This function performs a random walk in the hyperbolic upper half-plane
    % n_steps: Number of steps in the random walk
    % step_size: Fixed hyperbolic distance to move at each step

    lambda = 1.1;
    eps = 0;

    % Initialize the starting point
    z = 1i;  % Start at (0, 1) to avoid the real axis
    
    % Pre-allocate array for plotting
    points = zeros(n_steps + 1, 1);
    points(1) = z;
    
    % Iterative Process
    for k = 2:(n_steps + 1)
        % Get current point and the possible points to step to in a circle
        % around it
        z = points(k - 1);
        ranges = [0, 2*pi];
    
        % Determine the directions we can go without violating
        % quasi-geodesic lower bound condition
        for t_i = 1:min(k - 2, ceil(k - lambda*eps) - 1) % Note: t_i > lambda*eps - k <=> (k - t_i)/lambda - eps <= 0
            % Get the i-th previous segment and determine if 
            segment_i = GeodesicSegment(points(t_i), points(t_i + 1));
    
            % Determine the s-neighborhood around the segment, which is a
            % neighborhood within which we cannot guarentee stepping into
            % will preserve the quasi-geodesic lower bound
            lowerBd = (1/lambda) * abs(k-t_i) * step_size - eps;
            % TO-DO: Clean up this mess later
            if t_i ~= k - 2
                fat = (1 - 1/lambda) * step_size;
                max_s = max(lowerBd + fat, 0); % Fatten lower bound (might not be necessary)
    
                % Fixes bug when one circle can get stuck inside the other
                if step_size == max_s
                    max_s = max_s - 0.01;
                end
                
                intersection_i = intersections_of_point_and_segment_ngbhs(z, segment_i, ...
                    step_size, max_s);
                inside = false;
            else
                intersection_i = intersections_of_circles(z, 0.1, points(t_i), lowerBd);
                inside = true;
            end

            k_and_t = [k, t_i]
            lower_bound = lowerBd
            inters = double(intersection_i)
            
            range_i = getRangeTn(z,intersection_i, inside);
            range = double(range_i)
            ranges = [ranges;range_i];
        end
        merged_range = MergeRange(ranges);
    
        % generate the new point in the specified range
        points(k) = generateTn(z,merged_range,step_size);
    end
    
    % Plot the path
    figure;
    plot(real(points), imag(points), 'o-');
    xlabel('Real part');
    ylabel('Imaginary part');
    title('Random Walk in the Hyperbolic Upper Half Plane');
    axis equal;

    % Verify quasi-geodesicity
    verify_quasigeodesic(points, lambda, eps, step_size)
end

% calculate the center of the great geodesic circle given two points on it
function geoCenter = getGeoCenter(s,t)
    x1 = real(s);
    y1 = imag(s);
    x2 = real(t);
    y2 = imag(t);
        
    geoCenter = (x1^2 + y1^2 - x2^2 - y2^2) / (2 * (x1 - x2));
end

% calculate the slope of the tangent of a circle
function tanSlope = getTanSlope(point,center)

    u = point;
    t = center;
    c = getGeoCenter(u,t); % center of the great geodesic circle

    % Evaluate the derivative of circle with center c and point t
    % at the center of step size circle t
    tanSlope = (real(t) - real(c)) / (imag(c) - imag(t));
end


% calculate the angle between intersecting points
function rangeTn = getRangeTn(t, intersection, t_inside_boundary)
    if isempty(intersection)
        rangeTn = [0, 2*pi];
        return
    end
    u1 = intersection(1); 
    u2 = intersection(2); % TO-DO: Case for 4 intersections
    
    % calculate the slope of the tangent of a circle
    slope_1 = getTanSlope(u1,t);
    slope_2 = getTanSlope(u2,t);
    
    % calculate the angle by setting vertical line as angle zero
    t_x = real(t);
    u1_x = real(u1);
    u2_x = real(u2);
    if u1_x < t_x
        theta_1 = pi/2 + atan(slope_1); % negative slope gives negative atan  
    else
        theta_1 = 3*pi/2 + atan(slope_1);
    end

    if u2_x < t_x
        theta_2 = pi/2 + atan(slope_2); % negative slope gives negative atan  
    else
        theta_2 = 3*pi/2 + atan(slope_2); 
    end            
    
    % Calculate the angle range outside of boundary
    angDiff = abs(theta_1 - theta_2);
    if ~t_inside_boundary
        if angDiff < pi % then start from larger-angled arm
            rangeTn_start = max(theta_1,theta_2);
            rangeTn = [rangeTn_start, rangeTn_start + 2*pi - angDiff];
        else
            rangeTn_start = min(theta_1,theta_2);
            rangeTn = [rangeTn_start, rangeTn_start + angDiff];
        end
    else 
        if angDiff < pi % then start from larger-angled arm
            rangeTn_start = min(theta_1,theta_2);
            rangeTn = [rangeTn_start, rangeTn_start + angDiff];
        else
            rangeTn_start = max(theta_1,theta_2);
            rangeTn = [rangeTn_start, rangeTn_start + 2*pi - angDiff];
        end
    end
end



% merge multiple ranges
function mergedRange = MergeRange(ranges)
    % Initialize the intersection range
    intersection_min = -inf;
    intersection_max = inf;
    
    % Find the intersection of all ranges
    for i = 1:size(ranges, 1)
        intersection_min = max(intersection_min, ranges(i, 1));
        intersection_max = min(intersection_max, ranges(i, 2));
    end
    
    % Check if the intersection is valid
    if intersection_min >= intersection_max
        double(ranges)
        error('The ranges do not overlap.');
    end

    mergedRange = [intersection_min,intersection_max];
end

% generate bounded t_n using randomization
function tn = generateTn(t,range,step_size) % ?might be able to generate "geodesic"
    
    % generate a random angle within the specified range
    phi = range(1) + (range(2) - range(1)) * rand();

    % calculate the coordinates of the random point
    % NOTE: angles in range place 0 vertically up, while angles in the call
    % to this function place 0 to the right like normal polar coordinates.
    % Thus, we add -pi/2 for conversion.
    tn = get_point_along_direction(t, double(phi + pi/2), step_size);
end