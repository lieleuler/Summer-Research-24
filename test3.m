
%g1 = GeodesicSegment(1i, 3.1i);
%g2 = GeodesicSegment(0.8 + 0.3i, 2 + 1i);

g2 = GeodesicSegment(0.0000 + 1.0000i, 0.0783 + 0.9426i);
g1 = GeodesicSegment(0.0000 + 1.0000i, 0.4486 + 1.4928i);

%a1 = dist_H(1i, 1.0635 + 0.73241168i)
%a2 = dist_H(3.1i, 1.5092 + 0.98805258i)

%a1_with_dist_from_point = g2.dist_from_point(1i)
%a2_with_dist_from_point = g2.dist_from_point(3.1i)


tic
get_geodesic_directed_hausdorff(g1, g2)
toc



