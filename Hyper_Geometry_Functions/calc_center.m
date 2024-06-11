function c = calc_center(this)
    x1 = real(this.start_point);
    y1 = imag(this.start_point);
    x2 = real(this.end_point);
    y2 = imag(this.end_point);

    c = (x1^2 + y1^2 - x2^2 - y2^2) / (2 * (x1 - x2));
end
