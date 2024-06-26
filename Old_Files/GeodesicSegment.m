
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
            my_obj.start_point = p1;
            my_obj.end_point = p2;

            my_obj.is_line = (real(p1) == real(p2));
        end

        % == Geometry Methods == %
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
            y1 = min(imag(e1), imag(e2));
            y2 = max(imag(e1), imag(e2));
            transformed_point = (a*point + b) / (c*point + d);

            optimal_y = abs(transformed_point);

            % Hacky stuff to make it work with sym
            logic_intermediary = (y1 - optimal_y) * (y2 - optimal_y);
            % -1 if y1 <= optimal_y <= y2, 1 else
            logic_factor = abs(logic_intermediary) / logic_intermediary;
            dist = (1 - 1)/-2 * dist_H(transformed_point, optimal_y*1i) + ... 
                   (1 + 1)/2 * max(dist_H(transformed_point, e1), dist_H(transformed_point, e2));

            %if (y1 <= optimal_y && optimal_y <= y2) || (y2 <= optimal_y && optimal_y <= y1)
                %dist = dist_H(transformed_point, optimal_y*1i);
            %else
                %dist = max(dist_H(transformed_point, e1), dist_H(transformed_point, e2));
            %end
        end
        function point = travel_from_start(this, dist)
            [a, b, c, d] = this.find_flt_to_imag_axis();
            transformed_start_point = (a*this.start_point + b) / (c*this.start_point + d);
            point_dist_away = transformed_start_point * exp(dist);
            point = (d*point_dist_away - b) / (-c*point_dist_away + a);
        end
        % == Display Methods == %
        function plot(this, steps)
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
            plot(x_values, y_values)
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