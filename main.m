
DELTA = log(1 + sqrt(2));

iteration_num = 1;

lambda = 6;
eps = 1;
step_size = 0.1;
min_segment_splits = 1;

k = 0;

for i = 1:iteration_num
    tic
    [points, ranges, phi] = random_walk_hyperbolic(1000, lambda, eps, step_size, min_segment_splits);
    toc
    %if ~verify_quasigeodesic(points, lambda, eps, step_size, 20)
        %ranges
        %phi
        %k = k + 1;
    %end

    tic
    quasi_geodesic = points_to_quasigeodesic(points);
    true_geodesic = GeodesicSegment(points(1), points(length(points)));
    verification = verify(quasi_geodesic, true_geodesic, DELTA, lambda, eps, 15);
    toc
end

verification
%points
%bad_percent = k/iteration_num
