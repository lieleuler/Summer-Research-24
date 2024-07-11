
produces_quasi_geodesic_curves(6, 0.5, 0.1, 1, 1, 0);



function produces_quasi_geodesic_curves(lambda, eps, step_size, target_length, ...
    min_segment_splits, jump_chance)
    total_trials = 10;
    bad_quasi_geodesics = 0;
    for i = 1:total_trials
        points = random_walk_hyperbolic(lambda, eps, step_size, ...
                                        target_length, min_segment_splits, ...
                                        jump_chance);
        if ~verify_quasigeodesic(points, lambda, eps, step_size, 10)
            bad_quasi_geodesics = bad_quasi_geodesics + 1;
        end
    end
    bad_quasi_geodesics/total_trials
end