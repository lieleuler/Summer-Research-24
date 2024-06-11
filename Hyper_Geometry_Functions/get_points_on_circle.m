
function points = get_points_on_circle(center, radius, angle_divisions)
    angles = linspace(0, 2*pi, angle_divisions + 1);

    points = zeros(angle_divisions, 1);
    for t=1:angle_divisions
        points(t) = get_point_along_direction(center, angles(t), radius);
    end
end