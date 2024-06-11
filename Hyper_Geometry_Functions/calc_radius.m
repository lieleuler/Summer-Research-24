function r = calc_radius(this)
    x1 = real(this.start_point);
    y1 = imag(this.start_point);

    c = this.get_center_on_real_line();

    r = sqrt((x1 - c)^2 + y1^2);
end
