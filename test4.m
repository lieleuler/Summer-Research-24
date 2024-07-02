
tg6 = GeodesicSegment(0.2537 + 1.3207i,   0.3044 + 1.3702i);

qd15 = GeodesicSegment(0.3839 + 1.1867i, 0.3101 + 1.2859i);
qd16 = GeodesicSegment(0.3101 + 1.2859i, 0.4302 + 1.3389i);


z = tg6.travel_from_start_by_percent(0.5);

qd15.dist_from_point(z)
qd16.dist_from_point(z)