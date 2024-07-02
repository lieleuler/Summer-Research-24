
qg1 = GeodesicSegment(1i, 0.0783 + 0.9426i);
qg2 = GeodesicSegment(0.0783 + 0.9426i, 0.0552 + 1.0388i);
qg3 = GeodesicSegment(0.0552 + 1.0388i, 0.1445 + 1.0975i);

tg1 = GeodesicSegment(0.0195 + 1.0303i, 0.4486 + 1.4928i);
tg2 = GeodesicSegment(1i, 0.0195 + 1.0303i);

z = 0.0195 + 1.0303i;
w = 0.0783 + 0.9426i;


qg1.dist_from_point(z)
dist_H(z, w)

get_geodesic_directed_hausdorff(tg1, qg1)




