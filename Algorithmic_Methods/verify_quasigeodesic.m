
% Verifies if a collection of points generated by our algorithm is indeed a
% lambda-epsilon quasigeodesic

function is_quasigeodesic = verify_quasigeodesic(points, lambda, epsilon, ...
    step_size, subdivisions)
    num_points = size(points, 1);
    is_quasigeodesic = true;
    for i = 1:(num_points - 1)
        for j = i + 1:num_points
            p_i = points(i);
            p_j = points(j);

            t_i = step_size * (i - 1);
            t_j = step_size * (j - 1);

            if ~((t_j - t_i)/lambda - epsilon <= dist_H(p_i, p_j))
                i
                j
                is_quasigeodesic = false;
                return
            end

            segment_i = GeodesicSegment(points(i), points(i + 1));
            for i_sub = linspace(0, (subdivisions - 1) / subdivisions * step_size, subdivisions)
                p_i_sub = segment_i.travel_from_start(i_sub);
                for j_sub = linspace(0, (subdivisions - 1) / subdivisions * step_size, subdivisions)
                    if j < num_points
                        segment_j = GeodesicSegment(points(j), points(j+1));
                        p_j_sub = segment_j.travel_from_start(j_sub);
                    else
                        p_j_sub = p_j;
                        j_sub = 0;
                    end
                    if ~((t_j + j_sub - t_i - i_sub)/lambda - epsilon <= dist_H(p_i_sub, p_j_sub))
                        i + i_sub
                        j + j_sub
                        is_quasigeodesic = false;
                        return
                    end
                end
            end
        end
    end