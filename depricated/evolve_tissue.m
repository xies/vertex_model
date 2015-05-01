function dr = evolve_tissue(t,vert_coords,tis)
% Wrapper
t

eta = tis.parameters.viscosity;

tis = tis.evolve(vert_coords);
dr = tis.get_force;
dr = dr / eta;
dr = dr(:);

tis.draw('showContractile','showVectors',reshape(dr,numel(dr)/2,2));
drawnow
% tissueArray(counter) = Tissue;
% counter = counter + 1;

end