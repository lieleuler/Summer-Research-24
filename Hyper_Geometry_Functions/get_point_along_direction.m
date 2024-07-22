
function point = get_point_along_direction(start, phi, magnitude)
    % Define 0 degrees to be straight up in the upper half plane
    if mod(phi, 2*pi) == 0
        point = real(start) + exp(magnitude)*imag(start)*1i;
        return
    elseif mod(phi, pi) == 0
        point = real(start) + exp(-magnitude)*imag(start)*1i;
        return
    end

    % Perform hyperbolic rotation
    cot_phi_over_2 = cot(phi/2);
    exp_mag = exp(magnitude);
    point = real(start) + imag(start)*(1 + cot_phi_over_2*exp_mag*1i)/(cot_phi_over_2 - exp_mag*1i);
end