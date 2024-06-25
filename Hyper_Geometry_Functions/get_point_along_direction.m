
function point = get_point_along_direction(start, angle, magnitude)
    % Normalize angle to be in [0, 2pi)
    angle = mod(angle, 2*pi);

    if angle == pi/2
        point = real(start) + 1i*(imag(start) * exp(magnitude));
    elseif angle == 3*pi/2
        point = real(start) + 1i*(imag(start) * exp(-magnitude)); 
    else
        u = real(start);
        v = imag(start);
        center_on_real_axis = v/(tan(pi/2 - angle)) + u;
        mirrored_point = start + 2*(center_on_real_axis - u);
        geodesic = GeodesicSegment(start, mirrored_point);
        % Step along geodesic by the magnitude
        if xor(pi/2 <= angle && angle <= 3*pi/2, real(start) < real(mirrored_point))
            point = geodesic.travel_from_start(magnitude);
        else
            point = geodesic.travel_from_start(-magnitude);
        end
    end
end