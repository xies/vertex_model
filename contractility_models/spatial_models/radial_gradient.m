function C = radial_gradient(centroids,~,params)
%RADIAL_GRADIENT
% Sets a radial spatial gradient of contractility.
%
% USAGE: C = radial_gradient(centroids,t,params)
%
% INPUT: centroids - cell centroids
%        t - time (not used)
%        params(1) - Max contractility
%        params(2) - midpoint x
%        params(3) - midpoint y
%        params(4) - decay constant

A = params(1);
midline_x = params(2);
midline_y = params(3);
sigma = params(4);

cy = centroids(:,1); % YES, the first one-- matlab sucks.
cx = centroids(:,2);

r = sqrt((cy-midline_y).^2 + (cx - midline_x).^2);

C = A*exp(-(r).^2 ./ (2*sigma).^2 );