% calculate the directed hausdorff distance from quasi-geodesic to true geodesic
function quasi_geodesic_directed_hausdorff = get_quasi_geodesic_directed_hausdorff(quasi_geodesic,true_geodesic)
    quasi_geodesic_directed_hausdorff = 0; % initialize    
    for i = 1:length(quasi_geodesic)
        dist = get_geodesic_directed_hausdorff(quasi_geodesic(i),true_geodesic);
        if dist > quasi_geodesic_directed_hausdorff
            quasi_geodesic_directed_hausdorff = dist;
        end
    end
    assert(quasi_geodesic_directed_hausdorff > 0, "quasi_geodesic_directed_hausdorff is zero, something is wrong.")
end
