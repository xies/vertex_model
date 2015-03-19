function distances = compute_distances(verts, n)

y = verts(:,1); x = verts(:,2);
y_rep = repmat(y, 1, n);
x_rep = repmat(x, 1, n);
distances = sqrt((y_rep - y_rep').^2 + (x_rep - x_rep').^2);