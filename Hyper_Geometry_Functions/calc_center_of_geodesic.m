function c = calc_center_of_geodesic(g1, g2)
    x1 = real(g1);
    y1 = imag(g1);
    x2 = real(g2);
    y2 = imag(g2);

    c = (x1^2 + y1^2 - x2^2 - y2^2) / (2 * (x1 - x2));
end
