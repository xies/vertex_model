function varargout = run_model(init,params,contraction,OUT_DIR)
%RUN_MODEL
% Uses ODE23 (Runge-Kutta 2/3; Bockagi-Shampine) to solve a vertex model
% given by INITI, with parameters PARAM, and contraction setting
% CONTRACTION. Output is put into OUT_DIR.
%
% Usage:
%   run_model(initialization, parameters, contraction ...
%       OUT_DOR)
%
% INPUT: init - initialization structure
%           .initialize - lattice creator (e.g. create_hexagons)
%           .numX / Y - number of cells in each dimension
%           .hex_angle - angle of hexagon, (e.g. 'diagonal')
%           .t0 / tf - time span of solution
%           .cell_size - initial size of cells in microns^2
%           .abs_tol / rel_tol - tolerances for ODE23
%        parameters - structure of parameters to pass onto
%           Tissue.setParameteters
%        contraction
%           .ventral.box - ventral cell region
%           .ventral.alt_tension - alternative line tension for nonventral cells
%           .contract.model /parameter - function handle for how to set
%               contraction value (and parameters to pass to it)
%       OUT_DIR = output directory
%
% OUTPUT: optional: assembled tissue array
%
% See also: RUN_MODEL_FROM_CONFIG, Tissue/evolve, Tissue/step,
% assemble_model, ODE23

% Construct rough initial cells
initial_cells = init.initialize( init.model_params{:} );
[centroid_list, regions] = ...
    get_cents(initial_cells);
vertex_list = get_vertices(initial_cells);

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
        't1Threshold',params.t1Threshold ...
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
        't1Threshold',params.t1Threshold ...
        };
end

tis.setParameters(param_config{:});
T = toc;
display(['Parameter and connection matrices initialized in ' num2str(T) ' sec'])

% Set contractility
tic
% Activate "ventral fate"
cIDs = tis.getCellsWithinRegion(contraction.ventral.box);
tis.activateCell(cIDs, ...
    contraction.ventral.alt_tension*tis.parameters.areaElasticity);
% Set the value of contractility in each cell
tis.setContractilityModel( contraction );
T = toc;
display(['Contractility set in ' num2str(T) ' sec'])

% Add some jitter
tic
tis.jitterVertices(tis.parameters.jitterSize);
T = toc;
display(['Jitter added and contractility set in ' num2str(T) ' sec'])

% init.integration_method(tis,init,OUT_DIR);

% Need to make sure everything is gone in OUT_DIR that might
% conflict.
s = what(OUT_DIR);
if ~isempty(s.mat)
    error(['Please clear the contents of ' OUT_DIR]);
end

tisArr = tis.solve([init.t0 init.tf],OUT_DIR);
% Save initial tissue configuration for later assembling
if nargout > 0
    varargout{1} = tisArr;
end


end