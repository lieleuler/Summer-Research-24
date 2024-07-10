
classdef GeodesicSegment
    properties (Access = private)
        start_point
        end_point
        length
        is_line
        center
        radius
        angle
    end
    methods
        % == Constructor == %
        function my_obj = GeodesicSegment(p1, p2)
            if nargin > 0
                my_obj.start_point = p1;
                my_obj.end_point = p2;
    
                my_obj.is_line = (real(p1) == real(p2));
            end
        end

        % == Geometry Methods == %
        function midpoint = get_midpoint(this)
            midpoint = this.travel_from_start(this.get_length() / 2);
        end
        function p = travel_from_start_by_percent(this, percent)
            p = this.travel_from_start(this.get_length() * percent);
        end
        function new_geod = fractional_linear_transform(this, a, b, c, d) % mobius transformation
            new_start_point = (a * this.start_point + b) / (c * this.start_point + d);
            new_end_point = (a * this.end_point + b) / (c * this.end_point + d);
            new_geod = GeodesicSegment(new_start_point, new_end_point);
        end
        function [a, b, c, d] = find_flt_to_imag_axis(this) % find the coefs of the mobius transformation
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
        function dist = dist_from_point(this, point)
            [a, b, c, d] = this.find_flt_to_imag_axis();
            tranformed_geo = this.fractional_linear_transform(a, b, c, d);
            [e1, e2] = tranformed_geo.get_endpoints();
            y1 = imag(e1);
            y2 = imag(e2);
            transformed_point = (a*point + b) / (c*point + d);

            optimal_y = abs(transformed_point);

            if (y1 <= optimal_y && optimal_y <= y2) || (y2 <= optimal_y && optimal_y <= y1)
                dist = dist_H(transformed_point, optimal_y*1i);
            else
                dist = min(dist_H(transformed_point, e1), dist_H(transformed_point, e2));
            end
        end
        function point = travel_from_start(this, dist)
            [a, b, c, d] = this.find_flt_to_imag_axis();
            transformed_start_point = (a*this.start_point + b) / (c*this.start_point + d);
            transformed_end_point = (a*this.end_point + b) / (c*this.end_point + d);
            if imag(transformed_start_point) < imag(transformed_end_point)
                point_dist_away = transformed_start_point * exp(dist);
            else 
                point_dist_away = transformed_start_point * exp(-dist);
            end
            point = (d*point_dist_away - b) / (-c*point_dist_away + a);
        end
        function angle = get_angle_with_vertical(this, dist_along)
            if this.is_line
                angle = 0;
                return
            end

            c = this.get_center_on_real_line();
            p = this.travel_from_start(dist_along);
            p_x = real(p);
            p_y = imag(p);
            slope = -(p_x - c)/p_y;
            if xor(real(this.start_point) < real(this.end_point), dist_along < 0)
                % Case of traveling left
                angle = 3/2*pi + atan(slope);
            else
                % Case of traveling right
                angle = 1/2*pi + atan(slope);
            end
        end
        function does_intersect = intersects_geodesic(this, geod, can_extend_self, can_extend_geod)
            try
                does_intersect = 0 ~=...
                   height(this.intersections_with_geodesic(geod, ...
                   can_extend_self, can_extend_geod));
            catch ME
                if ME.identifier == "Geodesic:SameGeodesicError"
                    does_intersect = true;
                else
                    rethrow(ME)
                end
            end
        end
        function points = intersections_with_geodesic(this, geod, can_extend_self, can_extend_geod)
            points = [];

            this_min_x = min(real(this.start_point), real(this.end_point));
            this_max_x = max(real(this.start_point), real(this.end_point));
            [e1, e2] = geod.get_endpoints();
            geod_min_x = min(real(e1), real(e2));
            geod_max_x = max(real(e1), real(e2));

            c_1 = this.get_center_on_real_line();
            r_1 = this.get_radius_from_center();
            c_2 = geod.get_center_on_real_line();
            r_2 = geod.get_radius_from_center();

            if this.is_line
                if geod.is_line
                    if this_max_x == geod_max_x
                        error("Geodesic:SameGeodesicError", "The provided geodesics are identical")
                    end
                    return
                end
                if abs(r_2) > abs(c_2)
                    points = [sqrt(r_2^2 - c_2^2)*1i];
                end
                return
            end
            if geod.is_line
                if abs(r_1) > abs(c_1)
                    points = [sqrt(r_1^2 - c_1^2)*1i];
                end
                return
            end

            % If centers are equal, then geodesic can only intersect if
            % their radii are equal
            if c_1 == c_2
                if r_1 == r_2 && ...
                   (can_extend_self || can_extend_geod || ...
                   (this_max_x <= geod_min_x))
                    error("Geodesic:SameGeodesicError", "The provided geodesics are identical")
                end
                return
            end

            % Solve for intersection point across the whole geodesic, then
            % check if that point lies on the segments of each geodesic
            x = ((r_1^2 - r_2^2)/(c_2 - c_1) + c_1 + c_2)/2;
            
            % If x is out of bounds for both geodesics, this is sufficient
            % to say they can't intersect, since traveling along a geodesic 
            % either monotonically increases/decreases the x coordinate
            if ~ ( ((this_min_x <= x && x <= this_max_x) || can_extend_self) && ... 
                    ((geod_min_x <= x && x <= geod_max_x) || can_extend_geod) )
                return
            end

            % The x coord is already taken care of above. Only check needed
            % now is to see if the intersection is real (i.e. y is not
            % imaginary)
            y_squared = r_1^2 - (x - c_1)^2;
            if y_squared <= 0
                return
            else
                points = [points; x + 1i*sqrt(y_squared)];
            end
        end
        function R = get_region_of_point(this, p)
            [a, b, c, d] = this.find_flt_to_imag_axis();
            transformed_g = this.fractional_linear_transform(a, b, c, d);
            [e1, e2] = transformed_g.get_endpoints();
            transformed_p = (a*p + b) / (c*p + d);

            p_x = real(transformed_p);
            p_y = imag(transformed_p);
            e1_x = real(e1);
            e1_y = imag(e1);
            e2_x = real(e2);
            e2_y = imag(e2);
            
            e1_under = e1_y < e2_y;
            e2_under = ~e1_under;

            R = 2;
            if (e1_under && (p_x - e1_x)^2 + p_y^2 <= e1_y^2) || ...
               (~e1_under && (p_x - e1_x)^2 + p_y^2 >= e1_y^2)
                    R = 1;
            elseif (e2_under && (p_x - e2_x)^2 + p_y^2 <= e2_y^2) || ...
               (~e2_under && (p_x - e2_x)^2 + p_y^2 >= e2_y^2)
                    R = 3;
            end
        end
        % == Display Methods == %
        function plot(this, steps, color)
            if this.is_line
                x = real(this.start_point);
                y1 = imag(this.start_point);
                y2 = imag(this.end_point);
                x_values = linspace(x, x, steps);
                y_values = linspace(y1, y2, steps);
            else
                x1 = real(this.start_point);
                x2 = real(this.end_point);
                y1 = imag(this.start_point);
                y2 = imag(this.end_point);

                c = this.get_center_on_real_line();
                r = this.get_radius_from_center();

                theta1 = atan(y1/(x1 - c));
                if theta1 < 0
                    theta1 = pi + theta1;
                elseif theta1 == 0 && x1 < c
                    theta1 = pi;
                end
                theta2 = atan(y2/(x2 - c));
                if theta2 < 0
                    theta2 = pi + theta2;
                elseif theta2 == 0 && x2 < c
                    theta2 = pi;
                end
            
                t_values = linspace(0, 1, steps);
                test = r*t_values;
                x_values = r*cos((theta2-theta1)*t_values + theta1) + c;
                y_values = r*sin((theta2-theta1)*t_values + theta1);
            end
            plot(x_values, y_values, color)
        end
        % == Getters == %
        function [p1, p2] = get_endpoints(this)
            p1 = this.start_point;
            p2 = this.end_point;
        end
        function length = get_length(this)
            if not(exist(this.length, "var"))
                this.length = this.calc_length();
            end
            length = this.length;
        end
        function center = get_center_on_real_line(this)
            if not(exist(this.center, "var"))
                this.center = this.calc_center();
            end
            center = this.center;
        end
        function radius = get_radius_from_center(this)
            if not(exist(this.radius, "var"))
                this.radius = this.calc_radius();
            end
            radius = this.radius;
        end
    end
        % == Private Methods == %
    methods (Access = private, Hidden = true)
        function c = calc_center(this)
            x1 = real(this.start_point);
            y1 = imag(this.start_point);
            x2 = real(this.end_point);
            y2 = imag(this.end_point);
        
            c = (x1^2 + y1^2 - x2^2 - y2^2) / (2 * (x1 - x2));
        end
        function r = calc_radius(this)
            x1 = real(this.start_point);
            y1 = imag(this.start_point);

            c = this.get_center_on_real_line();
        
            r = sqrt((x1 - c)^2 + y1^2);
        end
        function l = calc_length(this)
            l = dist_H(this.start_point, this.end_point);
        end            
    end
end    