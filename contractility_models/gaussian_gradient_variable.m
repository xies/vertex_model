function C = gaussian_gradient_variable(centroids,params)

A = params(1);
midline = params(2);
sigma = params(3);
std = params(4);

cy = centroids(:,1); % YES, the first one-- matlab sucks.

C = A*exp(-(cy-midline).^2 ./ (2*sigma).^2 ) + ...
    std * randn(size(cy)) .* exp(-(cy-midline).^2 ./ (2*sigma).^2 );