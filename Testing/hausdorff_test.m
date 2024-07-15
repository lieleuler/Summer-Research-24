% test hausdorff distance
delta = log(1 + sqrt(2));
epsilon = 5;
lambda = 10;


lower = abs(1.7-0.56)/lambda - epsilon < dist_H(3+3.7i, 2.903+3.87i);
upper = dist_H(3+3.7i, 2.903+3.87i) < abs(1.7-0.56)*lambda + epsilon;
disp("lower cond: "+lower+"; upper cond: "+upper);


true_geodesic = GeodesicSegment(-6.9993 + 0.1i,6.9993 + 0.1i);
quasi_geodesic = [GeodesicSegment(-6.9993+ 0.1i,-1+6i),GeodesicSegment(-1+6i,-1+7i),...
    GeodesicSegment(-1+7i,3+3i),...
    GeodesicSegment(3+3i,3+4i),GeodesicSegment(3+4i,6.9993+ 0.1i)];
dist = verify(quasi_geodesic,true_geodesic,delta,lambda,epsilon,2);

disp("Hausdorff Distance: "+dist(2));

