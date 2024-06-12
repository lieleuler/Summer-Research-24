
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
        center_on_real_axis = v/(tan(pi/2 - angle)) - u;
        geodesic = GeodesicSegment(start, start + 2*(center_on_real_axis - u));
        % Step along geodesic by the step_size
        if ~(pi/2 < double(angle) && double(angle) < 3*pi/2)
            point = geodesic.travel_from_start(-magnitude);
        else
            point = geodesic.travel_from_start(magnitude);
        end
    end
end