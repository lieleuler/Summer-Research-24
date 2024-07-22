% manual verification of quasi-geodesic condition (very slow)
function is_quasigeodesic = verify_quasigeodesic(start_points, end_points, ...
    lambda, epsilon, step_size, subdivisions)
    num_points = length(start_points);
    is_quasigeodesic = true;
    for i = 1:(num_points - 1)
        for j = i + 1:num_points
            p_i = start_points(i);
            p_j = start_points(j);

            t_i = step_size * (i - 1);
            t_j = step_size * (j - 1);

            if ~((t_j - t_i)/lambda - epsilon <= dist_H(p_i, p_j) && dist_H(p_i, p_j) <= lambda*(t_j - t_i) + epsilon)
                i
                j
                [(t_j - t_i)/lambda - epsilon, dist_H(p_i, p_j), lambda*(t_j - t_i) + epsilon]
                is_quasigeodesic = false;
                %return
            end

            segment_i = GeodesicSegment(start_points(i), end_points(i));
            for i_sub = linspace(0, (subdivisions - 1) / subdivisions * step_size, subdivisions)
                p_i_sub = segment_i.travel_from_start(i_sub);
                for j_sub = linspace(0, (subdivisions - 1) / subdivisions * step_size, subdivisions)
                    if j < num_points
                        segment_j = GeodesicSegment(start_points(j), end_points(j));
                        p_j_sub = segment_j.travel_from_start(j_sub);
                    else
                        p_j_sub = p_j;
                        j_sub = 0;
                    end
                    if ~((t_j + j_sub - t_i - i_sub)/lambda - epsilon <= dist_H(p_i_sub, p_j_sub) && dist_H(p_i_sub, p_j_sub) <= lambda*(t_j + j_sub - t_i - i_sub) + epsilon)
                        %points
                        [i + i_sub, t_i + i_sub]
                        [j + j_sub, t_j + j_sub]
                        [(t_j + j_sub - t_i - i_sub)/lambda - epsilon, dist_H(p_i_sub, p_j_sub), lambda*(t_j + j_sub - t_i - i_sub) + epsilon]
                        is_quasigeodesic = false;
                        %return
                    end
                end
            end
        end
    end
end
