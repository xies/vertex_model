function tisArr = solve(tis,tspan,OUT_DIR)
% Uses ODE23 to solve model from given initial condition for a
% given time span.
%
% USAGE: tisArr = tis_init.solve(tpsna,OUT_DIR);

verts = tis.vert_coords;
save([OUT_DIR '/model_t_0.mat'],'tis');

opt = odeset('OutputFcn',@solver_output); % Suppress plotting
[T,Y] = ode23(@(t,y) tis.step(t,y),tspan,verts,opt);

csvwrite([OUT_DIR '/times.csv'],T);
csvwrite([OUT_DIR '/vertices.csv'],Y);

tisArr = assemble_model(OUT_DIR);

end % solve