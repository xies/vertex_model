function C = uniform(centroids,t,params)
%UNIFORM
% Sets the same contractility throughout ventral tissue.
%
% USAGE: C = uniform(centroids,t,params)
%
% INPUT: centroids - cell centroids
%        t - time (not used)
%        params(1) - Constant contractility

A = params(1);

C = A*ones(1,size(centroids,1));