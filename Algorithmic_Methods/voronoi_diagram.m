% ChatGPT example code for finding the voronoi diagram with vertices of the polygon as sites. Adaptation needed.

% Define the vertices of the polygon (example)
polygon_vertices = [0, 0; 6, 0; 6, 6; 4, 4; 2, 6; 0, 6];

% Compute Voronoi diagram using the polygon vertices as sites
[vx, vy] = voronoi(polygon_vertices(:,1), polygon_vertices(:,2));

% Create a polyshape object for the polygon
concave_polygon = polyshape(polygon_vertices);

% Plot the polygon
figure;
plot(concave_polygon, 'FaceColor', 'none', 'EdgeColor', 'k');
hold on;

% Plot the Voronoi diagram
plot(vx, vy, 'b-');

% Clip Voronoi regions to the polygon
[v, c] = voronoin(polygon_vertices);
for i = 1:length(c)
    if all(c{i} ~= 1)   % Skip points at infinity
        cell = v(c{i}, :);
        cell_poly = polyshape(cell);   % Create Voronoi cell polygon
        clipped_cell = intersect(concave_polygon, cell_poly);  % Clip to polygon
        plot(clipped_cell, 'FaceColor', 'r', 'FaceAlpha', 0.3);  % Plot the clipped cell
    end
end

% Plot the polygon vertices
plot(polygon_vertices(:,1), polygon_vertices(:,2), 'ko', 'MarkerFaceColor', 'k');

hold off;
title('Voronoi Diagram Clipped to Polygon');
