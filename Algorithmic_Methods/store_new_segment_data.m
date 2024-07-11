function [sub_segment_points, ...
          sub_segment_points_transformed, ...
          sub_segment_points_abcd_values] = ...
                store_new_segment_data(t_n, new_segment, step_size, segment_splits, ...
                                       sub_segment_points, ...
                                       sub_segment_points_transformed, ...
                                       sub_segment_points_abcd_values)

    for i = 0:(segment_splits)
        sub_i = i/segment_splits;
    
        sub_segment_end = new_segment.travel_from_start(sub_i*step_size);
        sub_segment_points(t_n, i + 1) = sub_segment_end;
    
        if i > 0
            sub_segment_start = sub_segment_points(t_n, i);
            sub_segment = GeodesicSegment(sub_segment_start, sub_segment_end);
            [a, b, c, d] = sub_segment.find_flt_to_imag_axis();
            [e1, e2] = sub_segment.fractional_linear_transform(a, b, c, d).get_endpoints();
            sub_segment_points_transformed(t_n, 2*i - 1) = e1;
            sub_segment_points_transformed(t_n, 2*i) = e2;
    
            sub_segment_points_abcd_values(t_n, 4*(i-1) + 1) = a;
            sub_segment_points_abcd_values(t_n, 4*(i-1) + 2) = b;
            sub_segment_points_abcd_values(t_n, 4*(i-1) + 3) = c;
            sub_segment_points_abcd_values(t_n, 4*(i-1) + 4) = d;
        end
    end
end
