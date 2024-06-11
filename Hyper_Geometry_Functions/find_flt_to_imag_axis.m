% find the coefs of the mobius transformation
function [a, b, c, d] = find_flt_to_imag_axis(this) 
    if this.is_line
        a = 1;
        b = -real(this.start_point);
        c = 0;
        d = 1;
    else
        cen = this.get_center_on_real_line();
        rad = this.get_radius_from_center();
        a = 1;
        b = -(cen + rad);
        c = 1;
        d = -(cen - rad);
    end
end
