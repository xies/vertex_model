function [E,gradient] = calculate_energy(verts, p)
% Calculate the current energy state. For use in minimizing energy.
% USAGE: [Energy,grad] = calculate_energy(verts, p_struct)
%
% INPUT: verts - Nv x 2 array of vertex coordinates
%        p - parameter structure
%         .fixed_cells - Nc x 1 logical array of which cells are fixed
%         .not_fixed_cells - Nc x 1 logical array of which cells are not fixed
%         .numV = Nv
%         .preferred_distances - "resting" length of all edges
%         .spring_constants - edge elasticity
%         .initial_verts - the initial vertices

verts = reshape(verts, length(verts)/2, 2); 
% put the fixed cells back in the verts array for distances
verts2 = p.initial_verts;
verts2(p.not_fixed_cells, :) = verts;
verts = verts2;
n = p.numV;

E = 0;

% calculate the distances between vertices
% tic,distances1 = compute_distances(verts, n); toc;
distances = squareform(pdist(verts));
distances(~p.connectivity) = NaN;
% diagonal elements = NaN
distances(logical(eye(size(distances)))) = NaN; 

% fix preferred_distances
preferred_distances = p.preferred_distances;
preferred_distances(~p.connectivity) = NaN;
preferred_distances(logical(eye(size(preferred_distances)))) = NaN; 

spring_constants = p.spring_constants;

% compute the energy
% IMPORTANT LINE
distance_fluct = spring_constants .* (distances - preferred_distances).^2;   
% % get rid of NaN elements for summing
% distance_fluct(isnan(distance_fluct)) = 0; 
% only keep half of the matrix (throws away diagonal as well)
distance_fluct = triu(distance_fluct); 
% distance_fluct = p.weights .* distance_fluct;
DISTANCES_TERM = nansum(distance_fluct(:));

E = E + DISTANCES_TERM;

% compute the gradient
gradient = compute_gradient(distances, preferred_distances, spring_constants, verts);

% remove the fixed cells from the gradient
gradient(p.fixed_cells, :) = [];

gradient = [gradient(:,1); gradient(:,2)];

end

% compute the Hessian

% hessian = compute_hessian(distances, preferred_distances, p.fixed_cells, verts, n);

% chol(hessian);

