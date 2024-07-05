
DELTA = log(1 + sqrt(2));

iteration_num = 1;

target_length = 10;
lambda = 5;
eps = 1;
step_size = 0.1;
min_segment_splits = 1;

k = 0;

for i = 1:iteration_num
    tic
    points = random_walk_hyperbolic(lambda, eps, step_size, target_length, min_segment_splits);
    toc
    num_points = length(points)
    %if ~verify_quasigeodesic(points, lambda, eps, step_size, 20)
        %"NOT VERIFIED"
    %end

    tic
    quasi_geodesic = points_to_quasigeodesic(points);
    true_geodesic = GeodesicSegment(points(1), points(length(points)));
    verification = verify(quasi_geodesic, true_geodesic, DELTA, lambda, eps, ceil(num_points/3));
    toc
end

verification
%points
%bad_percent = k/iteration_num
