function [E gradient hessian] = my_calc_energy(verts, p)
% change shape back to normal
verts = reshape(verts, length(verts)/2, 2); 

%% put the fixed cells back in the verts array for distances
verts2 = p.initial_verts;
verts2(p.not_fixed_cells, :) = verts;
verts = verts2;
n = p.numV;

E = 0;

% now cg is a temporary CellGraph with the current vertices
% cg.updateVertexCoords(verts);
% figure(1), imshow(cg.draw);

%% calculate the distances between vertices
distances = compute_distances(verts, n); 
distances(~p.connectivity) = NaN;
% diagonal elements = NaN
distances(logical(eye(size(distances)))) = NaN; 

%% fix preferred_distances
preferred_distances = p.preferred_distances;
preferred_distances(~p.connectivity) = NaN;
preferred_distances(logical(eye(size(preferred_distances)))) = NaN; 

spring_constants = p.spring_constants;

%% compute the energy
% IMPORTANT LINE
distance_fluct = spring_constants .* (distances - preferred_distances).^2;   
% get rid of NaN elements for summing
distance_fluct(isnan(distance_fluct)) = 0; 
% only keep half of the matrix (throws away diagonal as well)
distance_fluct = triu(distance_fluct); 
% distance_fluct = p.weights .* distance_fluct;
DISTANCES_TERM = sum(distance_fluct(:));

E = E + DISTANCES_TERM;

%% compute the gradient
gradient = compute_gradient(distances, preferred_distances, spring_constants, verts);

% remove the fixed cells from the gradient
gradient(p.fixed_cells, :) = [];

gradient = [gradient(:,1); gradient(:,2)];

%% compute the Hessian

% hessian = compute_hessian(distances, preferred_distances, p.fixed_cells, verts, n);

% chol(hessian);

