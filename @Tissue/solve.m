function tisArr = solve(tis,tspan,OUT_DIR)
% Uses ODE23 to solve model from given initial condition for a
% given time span.
%
% USAGE: tisArr = tis_init.solve(tpsna,OUT_DIR);

verts = tis.vert_coords;
tis.dir = OUT_DIR;
save([OUT_DIR '/model_t_0.mat'],'tis');

% Clear out the current contents of time.csv and vertices.csv
dlmwrite([OUT_DIR '/times.csv'],[]);
dlmwrite([OUT_DIR '/vertices.csv'],[]);

% Use custom output function
opt = odeset( ...
    'OutputFcn', @(t,y,flag) solver_output(t,y,flag,tis), ...
    'InitialStep', 1, ...
    'RelTol', 1e-4, ...
    'AbsTol', 1e-7 ...
    );

ode23(@(t,y) tis.step(t,y),tspan,verts,opt);

% csvwrite([OUT_DIR '/times.csv'],T);
% csvwrite([OUT_DIR '/vertices.csv'],Y);

tisArr = assemble_model(OUT_DIR);

end % solve