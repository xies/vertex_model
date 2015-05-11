function bogacki_shampine(tis,init,OUT_DIR)
% Performs Bogacki-Shampine integration on TISSUE model. Saves every step
% of the solution to OUT_DIR. Stopping criterion is both maxIter and an
% energy tolerance (relative)
%
% USAGE: bogacki_shampine(tis,init,OUT_DIR)

% Grab some initial parameters
STEPS = init.steps;
h = init.h0;
abs_tol = init.abs_tol;
rel_tol = init.rel_tol;

% Euler integration
E = 0;

for i = 1:STEPS
    
    % Copy previous tissue config
    tis_prev = tis; E_prev = E;
    
    tic
    if i == 1 % Only have to do this once, since first same as last
        % ---- K1 -----
        verts = tis_prev.vert_coords;
        k1 = tis_prev.get_force / tis_prev.parameters.viscosity;
    end
    
    % ---- K2 -----
    tmpTis = tis_prev.evolve( verts + k1 * h /2 );
    k2 = tmpTis.get_force / tmpTis.parameters.viscosity;
    
    % ---- K3 -----
    tmpTis = tis_prev.evolve( verts + k2 * h * 3/4 );
    k3 = tmpTis.get_force / tmpTis.parameters.viscosity;
    
    % Next step calculation
    displacements = ...
        ( 2/9 * k1 + 1/3 * k2 + 4/9 * k3) * h;
    tis = tis_prev.evolve( verts + displacements );
    
    % k1 for next step (used also for error estimation)
    k1 = tis.get_force / tmpTis.parameters.viscosity;
    
    Z = verts + h * (7/24 * k1 + 1/4 * k2 + 1/3 * k3 + 1/8 * k4);
    
    % Estimate next step size with error control
    scale = abs_tol + max( abs(tis.vert_coords),abs(Z) ) * rel_tol;
    error = max( abs(tis.vert_coords - Z),eps) / scale;
    h = h * error ^ (1/3);
    
    E = tis.get_energy;
    
    % Check if change in energy at this step is good enough
    if i > 1 && (abs(E - E_prev) < init.rel_tolerance * E(i) ...
            || abs(E - E_prev) < init.abs_tolerance )
        display(['Change in energy is ' num2str(E(i) - E(i-1))])
        break
    end
    
    % If change in energy is positive, then we're in an unstable situation
    % -> try to up-sample.
    if i>1 && E - E_prev > 0
        error('Unstable');
    end
    
    % Save current step in a .mat file
    T = toc;
    display([num2str(i) '-th time step (' num2str(T) ' sec)'])
    if nargin > 3
        SAVE_DIR = [OUT_DIR '_step_' num2str(i) '.mat'];
        save(SAVE_DIR,'tis');
        display(['Saved to: ' SAVE_DIR]);
    end
    
end