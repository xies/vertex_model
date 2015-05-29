function C = radial_gradient(centroids,params)
%RADIAL_GRADIENT

A = params(1);
midline_x = params(2);
midline_y = params(3);
sigma = params(4);

cy = centroids(:,1); % YES, the first one-- matlab sucks.
cx = centroids(:,2); % YES, the first one-- matlab sucks.

r = sqrt((cy-midline_y).^2 + (cx - midline_x).^2);

C = A*exp(-(r).^2 ./ (2*sigma).^2 );