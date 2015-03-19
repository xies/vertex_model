function gradient = compute_gradient(distances, preferred_distances, spring_constants, verts)

coeff = spring_constants.* (2*(distances - preferred_distances)./distances);
coeff(isnan(coeff)) = 0;

GRADIENT_DISTANCES_TERM = coeff * ones(size(verts)) .* ...
    verts - coeff * verts;

gradient = GRADIENT_DISTANCES_TERM;