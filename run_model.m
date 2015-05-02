function [tisArr,F] = run_model(init,params,contraction)
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
tisArr(1:STEPS+1) = Tissue;
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
tisArr(1) = tis;
E = zeros(1,STEPS);
for i = 1:STEPS
    
    tic
    verts = tis.vert_coords;
    displacements = tis.get_force ...
        /tis.parameters.viscosity * tis.parameters.lengthScale ...
        * tis.parameters.stepSize;
    verts = verts + displacements;
    
    tis = tis.evolve( verts );
    tisArr( i + 1 ) = tis;
    E(i) = tis.get_energy;
    
    if i>1 && abs(E(i) - E(i-1)) < init.tolerance
        display(['Change in energy is ' num2str(E(i) - E(i-1))])
        break
    end
    if i>1 && E(i) - E(i-1) > 0, error('Unstable regime!'); end
    
    T = toc;
    display([num2str(i) '-th time step (' num2str(T) ' sec)'])
    
end

tisArr(i+1:end) = [];
F = tisArr.movie('showContractility');

end