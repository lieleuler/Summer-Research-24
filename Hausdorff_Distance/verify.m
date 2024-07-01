% verify thm 17.1
function verification = verify(quasi_geodesic,true_geodesic,delta,lambda,epsilon)
    % calculate the theoretical upper bound
    R = get_R_value(delta,lambda,epsilon);
    
    % calculate the Hausdorff distance
    Hausdorff_distance = max(get_quasi_geodesic_directed_hausdorff(quasi_geodesic,true_geodesic),get_geodesic_directed_hausdorff(true_geodesic,quasi_geodesic));

    % verification
    verification = [R, Hausdorff_distance];
end
