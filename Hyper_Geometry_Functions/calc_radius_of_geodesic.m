function r = calc_radius_of_geodesic(g1, c)
    r = sqrt((real(g1) - c)^2 + imag(g1)^2);
end
