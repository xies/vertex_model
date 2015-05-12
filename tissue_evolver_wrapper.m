function dy = tissue_evolver_wrapper(t,y,tis,OUT_DIR)

y = reshape(y(:),[numel(y)/2 2]);

tis = tis.evolve(y,t);
dy = tis.get_force / tis.parameters.viscosity;
dy = dy(:);
save([OUT_DIR '_t_' num2str(t) '.mat'],'tis');

end