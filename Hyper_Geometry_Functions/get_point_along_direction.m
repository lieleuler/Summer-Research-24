
function point = get_point_along_direction(start, angle, magnitude)
    u = real(start);
    v = imag(start);
    center_on_real_axis = v*(pi/2 - tan(angle)) + u; 
    geodesic = GeodesicSegment(start, start + 2*(center_on_real_axis - u));
    % Step along geodesic by the step_size
    if 0 <= angle && angle < pi
        point = geodesic.travel_from_start(magnitude);
    else
        point = geodesic.travel_from_start(-magnitude);
    end
end
