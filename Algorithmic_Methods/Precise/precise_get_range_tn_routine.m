
function range_i = precise_get_range_tn_routine(z, e1, e2, step_size, raw_bound, fattened_bound, ...
    a, b, c, d, segment, is_lower_bound)
    
    [seg_start_pt, seg_end_pt] = segment.get_endpoints();
    seg_length = segment.get_length();
    % Case 1: Points are close enough to use regular algorithm
    if false && abs(log10(imag(seg_start_pt)/imag(z))) <= 5 && abs(log10(imag(seg_end_pt)/imag(z))) <= 5
        lb_intersections = precise_intersections_of_point_and_segment_ngbhs(z, ...
        e1, e2, step_size, fattened_bound, a, b, c, d);
        
        if ~isempty(lb_intersections)
            range_i = getRangeTn(z, lb_intersections, segment, step_size, ...
                                 fattened_bound, true);
        else
            range_i = [0, 2*pi];
        end
    % Case 2: Points are too far away from regular algorithm to be
    % precise... use sampling, which is slow and not analytic, but takes
    % advantage of the methods which we can make super precise
    else
        angle_precision = 100;
        % Get Bound
        if is_lower_bound
            bound = raw_bound + step_size + seg_length/2;
        else
            bound = raw_bound;
        end
        % Get Closest and Furthest Angles
        closest_point = seg_start_pt;
        if dist_H(z, seg_end_pt) <= dist_H(z, seg_start_pt)
            closest_point = seg_end_pt;
        end
        closest_angle = GeodesicSegment(z, closest_point).get_angle_with_vertical(dist_H(z, closest_point));
        furthest_angle = mod(closest_angle + pi, 2*pi);


        change_angles = zeros(2,2);
        iterations = ceil(log2(angle_precision));

        closest_angle_state = dist_H(get_point_along_direction(z, closest_angle, step_size), closest_point);
        furthest_angle_state = dist_H(get_point_along_direction(z, furthest_angle, step_size), closest_point);
        if closest_angle_state == furthest_angle_state
            if xor(closest_angle_state, is_lower_bound)
                range_i = [0, 2*pi];
            else
                range_i = [0, 0];
            end
            return
        end

        angle1 = closest_angle;
        angle2 = furthest_angle;
        if angle2 < angle1
            angle2 = angle2 + 2*pi;
        end
        for i = iterations
            mid_angle = (angle1 + angle2)/2;
            midpoint = get_point_along_direction(z, mid_angle, step_size);
            midpoint_state = dist_H(midpoint, closest_point) <= bound;
            if midpoint_state == closest_angle_state
                angle1 = mid_angle;
            else
                angle2 = mid_angle;
            end
        end
        change_angles(1, :) = [mod(angle1, 2*pi), closest_angle_state];

        angle1 = furthest_angle;
        angle2 = closest_angle;
        if angle2 < angle1
            angle2 = angle2 + 2*pi;
        end
        for i = iterations
            mid_angle = (angle1 + angle2)/2;
            midpoint = get_point_along_direction(z, mid_angle, step_size);
            midpoint_state = dist_H(midpoint, closest_point) <= bound;
            if midpoint_state == closest_angle_state
                angle2 = mid_angle;
            else
                angle1 = mid_angle;
            end
        end
        change_angles(2, :) = [mod(angle2, 2*pi), closest_angle_state];

        if xor(is_lower_bound, change_angles(1, 2))
            theta_1 = change_angles(1, 1);
            theta_2 = change_angles(2, 1) - 2*pi/floor(log2(angle_precision));
            if theta_2 < theta_1
                theta_2 = theta_2 + 2*pi;
            end
            range_i = [theta_1, theta_2];
        else
            theta_1 = change_angles(2, 1);
            theta_2 = change_angles(1, 1) - 2*pi/floor(log2(angle_precision));
            if theta_2 < theta_1
                theta_2 = theta_2 + 2*pi;
            end
            range_i = [theta_1, theta_2];
        end
    end
end




%         angle_samples = 12;
%         angles = linspace(0, 2*pi, angle_samples + 1); % Sample 0 twice (as 0 and 2pi) for algorithmic purposes
%         if is_lower_bound
%             bound = raw_bound + step_size + seg_length/2;
%         else
%             bound = raw_bound;
%         end
%         change_angles = zeros(0,2);
%         prev_state = false;
%         for angle = angles
%             point_at_angle = get_point_along_direction(z, angle, step_size);
%             within_bound = dist_H(point_at_angle, seg_start_pt) <= bound ...
%                            || dist_H(point_at_angle, seg_end_pt) <= bound;
%             if angle ~= 0 && within_bound ~= prev_state
%                 change_angles = [change_angles; [angle, within_bound]];
%             end
%             prev_state = within_bound;
%         end
%         if height(change_angles) == 0
%             if xor(prev_state, is_lower_bound)
%                 range_i = [0, 2*pi];
%             else
%                 % prev_state
%                 % is_lower_bound
%                 % bound
%                 % seg_start_pt
%                 % seg_end_pt
%                 % z_real = real(z)
%                 % z_imag = imag(z)
%                 % error("Sampling found range entirely inside lower bound neighborhood" + ...
%                 %     "or entirely outside upper bound neighborhood!")
%                 range_i = [0, 0];
%             end
%         elseif height(change_angles) == 2
%             if xor(is_lower_bound, change_angles(1, 2))
%                 range_i = [change_angles(1, 1), change_angles(2, 1) - 2*pi/angle_samples];
%             else
%                 range_i = [change_angles(2, 1), 2*pi + change_angles(1, 1) - 2*pi/angle_samples];
%             end
%         else
%             error("Weird amount of change angles found during sampling (" + height(change_angles) + ...
%                 " change angles)!")
