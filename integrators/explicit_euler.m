function tisArr = explicit_euler(tis,init,OUT_DIR)
% Performs explicit Euler integration on TISSUE model. Saves every step of
% the solution to OUT_DIR. Stopping criterion is both maxIter and an energy
% tolerance (relative)
% 
% USAGE: explicit_euler(tis,init,OUT_DIR)

STEPS = init.steps;

% Euler integration
E = 0;

for i = 1:STEPS
    
    % Copy previous tissue config
    tis_prev = tis; E_prev = E;
    
    tic
    verts = tis_prev.vert_coords;
    displacements = tis_prev.get_force ...
        / tis_prev.parameters.viscosity * tis_prev.parameters.lengthScale ...
        * tis_prev.parameters.stepSize;
    
    tis = tis_prev.evolve( verts + displacements);
    E = tis.energy
    
    % Check if change in energy at this step is good enough
    if i > 1 && (abs(E - E_prev) < init.rel_tol * E ...
            || abs(E - E_prev) < init.abs_tol )
        display(['Change in energy is ' num2str(E - E_prev)])
        break
    end
    
    % If change in energy is positive, then we're in an unstable situation
    % -> try to up-sample.
    stepSize = tis.parameters.stepSize;
    while i > 1 && E - E_prev > 0
        % Implement up-sampling by using smaller integration steps
        stepSize = stepSize / 2;
        display(['Trying smaller step size = ' num2str(stepSize)]);
        displacements = tis_prev.get_force ...
            / tis_prev.parameters.viscosity * tis_prev.parameters.lengthScale ...
            * stepSize;
        tis = tis_prev.evolve( verts + displacements );
        tis.parameters.stepSize = stepSize;
        E = tis.energy;
%         error('Unstable');
    end
    
    % Save current step in a .mat file
    T = toc;
    display([num2str(i) '-th time step (' num2str(T) ' sec)'])
    if nargin > 2
        SAVE_DIR = [OUT_DIR '_step_' num2str(i) '.mat'];
        save(SAVE_DIR,'tis');
        display(['Saved to: ' SAVE_DIR]);
    end
    
    tisArr(i) = tis;
    
end