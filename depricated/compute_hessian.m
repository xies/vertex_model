function hessian = compute_hessian(distances, preferred_distances, fixed_cells, verts, n)

y = verts(:,1); x = verts(:,2);
y_rep = repmat(y, 1, n);
x_rep = repmat(x, 1, n);

% hessian = compute_hessian(distances, p, x_rep, y_rep, n);
% split the matrix into 4 quadrants (TL, TR, BL, BR),
% and TR = BL' by the symmetry of the matrix.

% quadrant TL: x coordinates only
TL = -2 * ((x_rep - x_rep').^2 .* preferred_distances./(distances.^3) - ...
    (distances - preferred_distances)./distances);
TL(logical(eye(n))) = -sum(TL, 2);


% quadrant BR: y-coordinates only
BR = -2 * ((y_rep - y_rep').^2 .* preferred_distances./(distances.^3) - ...
    (distances - preferred_distances)./distances);
BR(logical(eye(n))) = -sum(BR, 2);


% quadrant TR'
TR = -2 * (x_rep - x_rep') .* (y_rep - y_rep') ...
    .* preferred_distances./(distances.^3);
TR(logical(eye(n))) = -sum(TR, 2);


% quadrant BL = TR'
BL = TR';

% remove the fixed cells from the Hessian
TL(fixed_cells, :) = []; TL(:, fixed_cells) = [];
BR(fixed_cells, :) = []; BR(:, fixed_cells) = [];
TR(fixed_cells, :) = []; TR(:, fixed_cells) = [];
BL(fixed_cells, :) = []; BL(:, fixed_cells) = [];

hessian = [TL TR; 
           BL BR];

hessian(isnan(hessian)) = 0;