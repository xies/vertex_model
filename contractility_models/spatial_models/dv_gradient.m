function C = dv_gradient(centroids,~,params)
%DV_GRADIENT
% Sets a spatial gradient of contractility in the 'DV' direction.
%
% USAGE: C = dv_gradient(centroids,t,params)
%
% INPUT: centroids - cell centroids
%        t - time (not used)
%        params(1) - Max contractility
%        params(2) - midpoint y
%        params(3) - decay constant

A = params(1);
midline = params(2);
sigma = params(3);

cy = centroids(:,1); % YES, the first one-- matlab sucks.

C = A*exp(-(cy-midline).^2 ./ (2*sigma).^2 );