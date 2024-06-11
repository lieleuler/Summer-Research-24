function point = travel_from_start(this, dist)
    [a, b, c, d] = this.find_flt_to_imag_axis();
    transformed_start_point = (a*this.start_point + b) / (c*this.start_point + d);
    point_dist_away = transformed_start_point * exp(dist);
    point = (d*point_dist_away - b) / (-c*point_dist_away + a);
end
