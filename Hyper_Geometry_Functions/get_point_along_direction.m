
function point = get_point_along_direction(start, phi, magnitude)
    % Define 0 degrees to be straight up in the upper half plane
    point_at_0_degrees = real(start) + (imag(start) * exp(magnitude))*1i;
    if mod(phi, 2*pi) == 0
        point = point_at_0_degrees;
        return
    end

    % Perform hyperbolic rotation
    cot_phi_over_2 = cot(phi/2);
    point = ((real(start) - imag(start)*cot_phi_over_2)*point_at_0_degrees - abs(start)^2)...
            /...
            (point_at_0_degrees - ((real(start) + imag(start)*cot_phi_over_2)));
end