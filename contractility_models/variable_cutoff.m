function C = variable_cutoff(centroids,params)

A = params(1);
B = params(2);

C = A*ones(1,size(centroids,1)) + B*randn(1,size(centroids,1));