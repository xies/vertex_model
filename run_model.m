function run_model(init,params,contraction,OUT_DIR)
% Wrapper for running vertex model

% Construct rough initial cells
initial_cells = init.initialize( init.model_params{:} );
[centroid_list, regions] = ...
    get_cents(initial_cells);
vertex_list = get_vertices(initial_cells);

STEPS = init.steps;

% Initialize Tissue
tic
tis = Tissue(regions,vertex_list,centroid_list,params.connectivity);
T = toc;
display(['Tissue initialized in ' num2str(T) ' sec'])

% Set parameters
A0 = mean([tis.getCells.area]);
P0 = mean([tis.getCells.perimeter]);
l = P0/6; % lattice length_scale
um_per_px = sqrt(init.cell_size/A0); % pixel size

if params.dimensionless
    param_config = {...
        'dimensonless', true, ...
        'lengthScale', l, ...
        'forceScale', params.forceScale, ...
        'lineTension', params.lineTension, ...
        'lineAnisotropy', 1, ...
        'areaElasticity', params.areaElasticity*l^3/params.forceScale, ...
        'perimElasticity', params.perimElasticity*l/params.forceScale, ...
        'targetArea',sqrt(3)*3/2, ...
        'viscosity', params.dragCoeff * params.forceScale / l, ...
        'connectivity', params.connectivity, ...
        'stepSize', params.timeStep, ...
        'jitterSize', params.jitterSize, ...
        'um_per_px', um_per_px, ...
        'dt_per_frame', 10, ...
        't1Threshold',3 ...
        };
else
    param_config = {...
        'dimensonless', false, ...
        'lengthScale', 1, ...
        'forceScale', params.forceScale, ...
        'lineTension', params.lineTension, ...
        'lineAnisotropy', 1, ...
        'areaElasticity', params.areaElasticity, ...
        'perimElasticity', params.perimElasticity, ...
        'viscosity', params.dragCoeff, ...
        'connectivity', params.connectivity, ...
        'stepSize', params.timeStep, ...
        'jitterSize', params.jitterSize, ...
        'um_per_px', um_per_px, ...
        'dt_per_frame', 10 ...
        't1Threshold',0 ...
        };
end

tis = tis.setParameters(param_config{:});
T = toc;
display(['Parameter and connection matrices initialized in ' num2str(T) ' sec'])

% Set contractility
tic
% Activate "ventral fate"
cIDs = tis.getCellsWithinRegion(contraction.ventral.box);
tis = tis.activateCell(cIDs, ...
    contraction.ventral.alt_tension*tis.parameters.areaElasticity);
% Set the value of contractility in each cell
tis = tis.setContractilityModel( ...
    contraction.contractility.model,contraction.contractility.params);
T = toc;
display(['Contractility set in ' num2str(T) ' sec'])

% Add some jitter
tic
tis = tis.jitterVertices(tis.parameters.jitterSize);
T = toc;
display(['Jitter added and contractility set in ' num2str(T) ' sec'])

% Euler integration
E = zeros(1,STEPS);
for i = 1:STEPS
    
    % Copy previous tissue config
    tis_prev = tis;
    
    tic
    verts = tis_prev.vert_coords;
    displacements = tis_prev.get_force ...
        / tis_prev.parameters.viscosity * tis_prev.parameters.lengthScale ...
        * tis_prev.parameters.stepSize;
    
    tis = tis_prev.evolve( verts + displacements);
    E(i) = tis.get_energy;
    
    % Check if change in energy at this step is good enough
    if i > 1 && abs(E(i) - E(i-1)) < init.tolerance * E(i)
        display(['Change in energy is ' num2str(E(i) - E(i-1))])
        break
    end
    
    % If change in energy is positive, then we're in an unstable situation
    % -> try to up-sample.
    stepSize = tis.parameters.stepSize;
    if i>1 && E(i) - E(i-1) > 0
        % Implement up-sampling by using smaller integration steps
        stepSize = stepSize / 2;
        displacements = tis_prev.get_force ...
            / tis_prev.parameters.viscosity * tis_prev.parameters.lengthScale ...
            * stepSize;
        tis = tis_prev.evolve( verts + displacements );
        tis.parameters.stepSize = stepSize;
        E(i) = tis.get_energy;
%         error('Unstable');
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

end