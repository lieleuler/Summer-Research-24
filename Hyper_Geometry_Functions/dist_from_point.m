function dist = dist_from_point(this, point)
    [a, b, c, d] = this.find_flt_to_imag_axis();
    tranformed_geo = this.fractional_linear_transform(a, b, c, d);
    [e1, e2] = tranformed_geo.get_endpoints();
    y1 = min(imag(e1), imag(e2));
    y2 = max(imag(e1), imag(e2));
    transformed_point = (a*point + b) / (c*point + d);

    optimal_y = abs(transformed_point);

    % Hacky stuff to make it work with sym
    logic_intermediary = (y1 - optimal_y) * (y2 - optimal_y);
    % -1 if y1 <= optimal_y <= y2, 1 else
    logic_factor = abs(logic_intermediary) / logic_intermediary;
    dist = (1 - 1)/-2 * dist_H(transformed_point, optimal_y*1i) + ... 
           (1 + 1)/2 * max(dist_H(transformed_point, e1), dist_H(transformed_point, e2));

    %if (y1 <= optimal_y && optimal_y <= y2) || (y2 <= optimal_y && optimal_y <= y1)
        %dist = dist_H(transformed_point, optimal_y*1i);
    %else
        %dist = max(dist_H(transformed_point, e1), dist_H(transformed_point, e2));
    %end
end
