% find the coefs of the mobius transformation

function [a, b, c, d] = find_flt_to_imag_axis(g1, g2) 
    if real(g1) == real(g2)
        a = 1;
        b = -real(g1);
        c = 0;
        d = 1;
    else
        cen = calc_center_of_geodesic(g1, g2);
        rad = calc_radius_of_geodesic(g1, cen);
        a = 1;
        b = -(cen + rad);
        c = 1;
        d = -(cen - rad);
    end
end
