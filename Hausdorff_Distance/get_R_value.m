% calculation for R in THM 1.7

function R_value = get_R_value(delta,lambda,epsilon)
    del = delta;
    lam = lambda;
    eps = epsilon;
    
    k_1 = lam * (lam + eps);
    k_2 = (2*lam * (lam + eps) + 3) * (lam + eps);
    
    a = (1/del) * log(2);
    b = (-6*k_1 - 2) * exp(log(2)/del);
    c = k_2 * exp(log(2)/del);
    
    % use the W_-1 branch of the Lambert-W function
    w = lambertw(-1, (a/b) * exp(a*c/b));
    
    D_0 = c/b - w/a;
    R_value = D_0 * (1 + k_1) + k_2 / 2;
end



% calculate the shortest distance between two geodesic segments
function dist_geod_to_geod = get_dist_geod_to_geod(geod_1, geod_2)
    % To be completeed
end

% calculate the Hausdorff distance between the geodesic and quasi-geodesics
function dist_Hausdorff = get_dist_Hausdorff(geod_segments, startpoint, endpoint)
    geodesic = GeodesicSegment(startpoint, endpoint);
distances = get_dist_geod_to_geod(geod_segments,geodesic); % a set of shortest paths from quasi-geod to geod
    dist_Hausdorff = max(distances);
end
