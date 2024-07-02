c = center(geodesic_2); % calculate the center of the circle on which geodesic_2 lies
    r = radius(geodesic_2); % calculate the radius of the circle on which geodesic_2 lies

    syms si % angle at center of d2 from x-axis to the chosen point p on s

    % calculate the range of si based on end points of the segement
    start_pt_1 = get_start_point(geodesic_1);
    end_pt_1 = get_end_point(geodesic_1);
    
    range_si = sort(atan(-start_pt_(2)/c),atan(-end_pt_1(2)/c));

    % calculate the range of theta based on end points of the segement
    start_pt_2 = get_start_point(geodesic_2);
    end_pt_2 = get_end_point(geodesic_2);

    range_theta = sort(atan((start_pt_2(1)-c)/start_pt_2(2)),atan((end_pt_2(1)-c)/end_pt_2(2)));

    % take derivative of distance function
    % Define the expression for x
    x = sqrt(-(c^2 + 2*c*r*cos(si))/(r^2));
    
    % Define the expression for g(si)
    g_si = 1 - (c^2 + r^2*x^2 + 2*c*r*x*(x*cos(si) + sqrt(1-x^2)*sin(si))) / (c*r*x*sqrt(1-x^2)*sin(si));
    
    % Compute the derivative of g(si) with respect to si
    g_prime = diff(g_si, si);
    
    % Convert the symbolic expression to a MATLAB function
    g_prime_func = matlabFunction(g_prime, 'Vars', si);

    % Use numerical solver to find si in the range range_si
    si_roots = vpasolve(g_prime_func, si, range_si);
    
    % values at si_roots
    x_si = sqrt(-(c^2 + 2*c*r*cos(si_roots))/(r^2)); % evaluate numerical value of x = cos(theta) at si_roots
    theta_si = acos(x_si); % evaluate numerical value of theta at si_roots
    
    % output MD
    if isempty(si_roots) % no critical point in the specified range_si
        minimal_distance = min(dist_H(start_pt_1,start_pt_2),dist_H(start_pt_1,end_pt_2),dist_H(start_pt_2,end_pt_1),dist_H(end_pt_1,end_pt_2));
    elseif theta_si > range_theta(2) || theta_si < range_theta(1) % theta_si outside the range_theta 
        minimal_distance = min(dist_H(start_pt_1,start_pt_2),dist_H(start_pt_1,end_pt_2),dist_H(start_pt_2,end_pt_1),dist_H(end_pt_1,end_pt_2));
    else % normal case 
        minimal_distance = acosh(1 - (c^2 + r^2*x^2 + 2*c*r*x*(x*cos(si_roots) + sqrt(1-x^2)*sin(si_roots))) / (c*r*x*sqrt(1-x^2)*sin(si_roots)));
    end
