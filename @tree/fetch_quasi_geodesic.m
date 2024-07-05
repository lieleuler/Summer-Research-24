% fetch the quasi-geodesic path
function quasi_geodesic_branch = fetch_quasi_geodesic(rwalk_tree, current_step, current_node_id)
    quasi_geodesic_branch = [];
    node_id = current_node_id;
    
    for t = current_step:1
        parent_id = getparent(node_id);
        quasi_geodesic_branch = [quasi_geodesic_branch, GeodesicSegment(rwalk_tree.get(node_id), rwalk_tree.get(parent_id))];
        node_id = parent_id;
    end
    
    quasi_geodesic_branch = quasi_geodesic_branch(end:-1:1); % reversing order
end