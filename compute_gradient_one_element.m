function gradient = compute_gradient_one_element(elem, distances, preferred_distances, spring_constants, verts)

distances = distances(:, elem);
preferred_distances = preferred_distances(:, elem);
spring_constants = spring_constants(:, elem);

coeff = spring_constants .* (2*(distances - preferred_distances)./distances);
coeff(isnan(coeff)) = 0;

GRADIENT_DISTANCES_TERM = repmat(coeff, 1, 2) .* ...
    (repmat(verts(elem, :), size(verts, 1), 1) - verts);

gradient = sum(GRADIENT_DISTANCES_TERM, 1);