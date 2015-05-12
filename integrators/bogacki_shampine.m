function tisArr = bogacki_shampine(tis,init,OUT_DIR)
% Performs Bogacki-Shampine integration on TISSUE model. Saves every step
% of the solution to OUT_DIR. Stopping criterion is both maxIter and an
% energy tolerance (relative)
%
% USAGE: bogacki_shampine(tis,init,OUT_DIR)

% Grab some initial parameters
STEPS = init.steps;
h = tis.parameters.stepSize;
abs_tol = init.abs_tol;
% rel_tol = init.rel_tol;

% Euler integration
E = tis.energy;

for i = 1:STEPS
    
    % Copy previous tissue config
    tis_prev = tis; E_prev = E;
    
    tic
    if i == 1 % Only have to do this once, since first same as last
        % ---- K1 -----
        verts = tis_prev.vert_coords;
        k1 = tis_prev.get_force / tis_prev.parameters.viscosity * h;
    else
        k1 = k4;
    end
    display(['Prev energy: ' num2str(tis_prev.energy - E_prev)]);
    
    % ---- K2 -----
    tmpTis = tis_prev.evolve( verts + k1 / 2 );
    k2 = tmpTis.get_force / tmpTis.parameters.viscosity * h;
    display(['K2 energy: ' num2str(tmpTis.energy - E_prev)]);
    
    % ---- K3 -----
    tmpTis = tis_prev.evolve( verts + k2 * 3/4 );
    k3 = tmpTis.get_force / tmpTis.parameters.viscosity * h;
    display(['K3 energy: ' num2str(tmpTis.energy - E_prev)]);
    
    % Next step calculation
    displacements = ...
        ( 2 * k1 + 3 * k2 + 4 * k3 ) / 9;
    tis = tis_prev.evolve( verts + displacements );
    
    % ---- k4 = k1 for next step (used also for error estimation)
    k4 = tis.get_force / tis.parameters.viscosity * h;
    
    % ---- err ----
    err = abs(5 * k1 - 6 * k2 - 8 * k3 + 9 * k4) / 72;
    
    % Estimate next step size with error control
    err = max( cat(1,err(:), eps) );
%     if err > abs_tol
%         h = h * max(abs_tol/err) ^ (1/2);
%     else
%         h = h * max(abs_tol/err) ^ (1/3);
%     end
%     h = .9 * h * max(max(abs_tol * tis.vert_coords)) / err;
    % Record next time step
    tis.parameters.stepSize = h;
    
    E = tis.energy
    display(['Updated energy: ' num2str(E)]);
    
    % Check if change in energy at this step is good enough
%     if i > 1 && (abs(E - E_prev) < init.rel_tol * E ...
%             || abs(E - E_prev) < init.abs_tol )
%         display(['Change in energy is ' num2str(E - E_prev)])
%         break
%     end
    
    % If change in energy is positive, then we're in an unstable situation
%     if i>1 && E - E_prev > 0, error('Unstable'); end
    
%     tis.draw('showContractile','showVectors',displacements*10);
%     drawnow;
    
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