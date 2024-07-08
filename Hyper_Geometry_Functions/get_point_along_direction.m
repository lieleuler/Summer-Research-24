
function point = get_point_along_direction(start, angle, magnitude)
    % Normalize angle to be in [0, 2pi)
    angle = mod(angle, 2*pi);

    if angle == 0
        point = real(start) + 1i*(imag(start) * exp(magnitude));
    elseif angle == pi
        point = real(start) + 1i*(imag(start) * exp(-magnitude)); 
    else
        u = real(start);
        v = imag(start);
        if angle ~= pi/2 && angle ~= 3*pi/2
            center_on_real_axis = -v/(tan(angle)) + u;
            % Mirror start across center of geodesic
            other_point = start + 2*(center_on_real_axis - u); 
        else
            other_point = u + v;
        end
        geodesic = GeodesicSegment(start, other_point);
        % Step along geodesic by the magnitude
        if xor(0 <= angle && angle <= pi, real(start) < real(other_point))
            point = geodesic.travel_from_start(magnitude);
        else
            point = geodesic.travel_from_start(-magnitude);
        end
    end
end