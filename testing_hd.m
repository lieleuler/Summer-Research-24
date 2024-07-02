
points = [
   0.0000 + 1.0000i
   0.0783 + 0.9426i
   0.0552 + 1.0388i
   0.1445 + 1.0975i
   0.1930 + 1.2016i
   0.2965 + 1.1462i
   0.4111 + 1.1460i
   0.5247 + 1.1354i
   0.4589 + 1.0483i
   0.5410 + 0.9880i
   0.4500 + 1.0320i
   0.5463 + 0.9997i
   0.5244 + 1.1024i
   0.4179 + 1.0788i
   0.3839 + 1.1867i
   0.3101 + 1.2859i
   0.4302 + 1.3389i
   0.5515 + 1.2886i
   0.4390 + 1.2319i
   0.4807 + 1.3542i
   0.4486 + 1.4928i
];

n = length(points);

true_geodesic = GeodesicSegment(points(1), points(n));

quasi_geodesic = [];
for i = 1:(n-1)
    quasi_geodesic = [quasi_geodesic, GeodesicSegment(points(i), points(i + 1))];
end

tic
[dist, intersections_segments, corresponding_segments] = get_polyline_directed_hausdorff_distance(true_geodesic, quasi_geodesic);
toc

dist = dist

%[e1, e2] = best_intersection.get_endpoints()

%plot(real(points), imag(points), "k")
%hold on
%GeodesicSegment(points(1), points(n)).plot(100, "k")
%hold on

colors = ["r", "g", "b", "c", "m", "y"];
colors_len = length(colors);

for i = 1:length(intersections_segments)
    color = colors(mod(i - 1, colors_len) + 1);
    intersections_segments(i).plot(100, color)
    hold on
    corresponding_segments(i).plot(100, color)

    [e1, e2] = intersections_segments(i).get_endpoints();
    [e3, e4] = corresponding_segments(i).get_endpoints();
    [i, e1, e2, e3, e4]
end

for i = 1:length(intersections_segments)
    is = intersections_segments(i);
    cs = corresponding_segments(i);

    closest_dist = Inf;
    for j = 1:length(quasi_geodesic)
        new_dist = get_minimal_distance(is, quasi_geodesic(j));
        if new_dist <= closest_dist
            closest_dist = new_dist;
            closest_seg = quasi_geodesic(j);
        end
    end

    i
    if closest_seg.get_endpoints() ~= cs.get_endpoints()
        "FATAL D: D: D:"
        [closest_seg.get_endpoints(); cs.get_endpoints()]
    end
end

% i = 6

% "FATAL D: D: D:"

%ans =
%
%   0.3101 + 1.2859i
%   0.3839 + 1.1867i

