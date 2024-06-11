function new_geod = fractional_linear_transform(this, a, b, c, d) % mobius transformation
    new_start_point = (a * this.start_point + b) / (c * this.start_point + d);
    new_end_point = (a * this.end_point + b) / (c * this.end_point + d);
    new_geod = GeodesicSegment(new_start_point, new_end_point);
end
