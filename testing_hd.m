
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

points = [
   0.0000 + 1.0000i
   0.0781 + 0.9423i
  -0.0074 + 0.9070i
   0.0301 + 0.9943i
   0.0760 + 1.0877i
   0.0085 + 1.0076i
  -0.0924 + 1.0126i
  -0.1750 + 0.9588i
  -0.2211 + 1.0479i
  -0.3081 + 0.9944i
  -0.3712 + 1.0764i
  -0.4481 + 1.1575i
  -0.5463 + 1.1017i
  -0.4387 + 1.0827i
  -0.5390 + 1.1293i
  -0.6522 + 1.1341i
  -0.5740 + 1.2222i
  -0.5542 + 1.1075i
  -0.6600 + 1.0796i
  -0.5823 + 1.0098i
  -0.6456 + 0.9360i
  -0.5582 + 0.9070i
  -0.5911 + 0.8268i
  -0.5267 + 0.8831i
  -0.4438 + 0.8565i
  -0.4945 + 0.7916i
  -0.5624 + 0.8364i
  -0.6428 + 0.8170i
  -0.7228 + 0.8036i
  -0.7874 + 0.8556i
  -0.7975 + 0.7747i
  -0.7688 + 0.7065i
  -0.7271 + 0.6529i
  -0.7350 + 0.7211i
  -0.6651 + 0.7426i
  -0.6889 + 0.8168i
  -0.7629 + 0.7858i
  -0.8398 + 0.8063i
  -0.9167 + 0.7857i
  -0.8572 + 0.8412i
  -0.8954 + 0.7703i
  -0.9146 + 0.6994i
  -0.8455 + 0.6913i
  -0.9064 + 0.7276i
  -0.8354 + 0.7475i
  -0.7725 + 0.7920i
  -0.6965 + 0.8185i
  -0.7701 + 0.7866i
  -0.8399 + 0.8270i
  -0.7653 + 0.8671i
  -0.8274 + 0.9322i
  -0.7885 + 1.0217i
  -0.7605 + 1.1253i
  -0.6881 + 1.0445i
  -0.5891 + 1.0160i
  -0.6117 + 0.9219i
  -0.5537 + 0.8546i
  -0.4698 + 0.8759i
  -0.4879 + 0.9661i
  -0.4770 + 1.0671i
  -0.4378 + 0.9730i
  -0.5298 + 0.9455i
  -0.4393 + 0.9220i
  -0.5081 + 0.8649i
  -0.4302 + 0.8312i
  -0.3509 + 0.8100i
  -0.4298 + 0.8332i
  -0.3473 + 0.8502i
  -0.2896 + 0.7918i
  -0.3688 + 0.8003i
  -0.2947 + 0.7738i
  -0.3694 + 0.7981i
  -0.3774 + 0.7225i
  -0.3069 + 0.7425i
  -0.3717 + 0.7827i
  -0.2946 + 0.7724i
  -0.2885 + 0.6991i
  -0.3566 + 0.6861i
  -0.3235 + 0.6293i
  -0.3194 + 0.5695i
  -0.2898 + 0.6211i
  -0.2577 + 0.5710i
  -0.3011 + 0.5366i
  -0.3042 + 0.5930i
  -0.2601 + 0.5561i
  -0.2610 + 0.6146i
  -0.3138 + 0.6493i
  -0.3642 + 0.6115i
  -0.3658 + 0.6758i
  -0.3957 + 0.6185i
  -0.3786 + 0.5620i
  -0.4102 + 0.5183i
  -0.4336 + 0.4745i
  -0.4736 + 0.4512i
  -0.4491 + 0.4914i
  -0.4119 + 0.4616i
  -0.4575 + 0.4566i
  -0.4448 + 0.5028i
  -0.3992 + 0.4840i
  -0.3735 + 0.4453i
  -0.4098 + 0.4735i
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
[e1, e2] = best_segment.get_endpoints()
[e3, e4] = best_intersection.get_endpoints()

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
    [i, e1, e2, e3, e4];
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

