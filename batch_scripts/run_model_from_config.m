function tisArr = run_model_from_config(CONFIG_DIR,OUT_DIR)

T = readtable(CONFIG_DIR,'ReadRowNames',true);
get_field = @(fieldname) cell2mat(T(fieldname,:).Value);

% Set parameters
initialization.initialize = @create_hexagons;
numX = str2double(get_field('num_cell_x'));
numY = str2double(get_field('num_cell_y'));
hex_angle = get_field('hex_angle');
initialization.model_params = {hex_angle,numX,numY};
initialization.t0 = str2double(get_field('t0')); % time span
initialization.tf = str2double(get_field('tf'));
initialization.abs_tol = str2double(get_field('abs_tol'));
initialization.rel_tol = str2double(get_field('rel_tol'));
l = 512 / max(numX,numY) / 2;
initialization.cell_size = ...
    str2double(get_field('cell_size')); % number of constriction steps

params.dimensionless = str2double(get_field('dimensionless'));
params.areaElasticity = str2double(get_field('kA'));
params.perimElasticity = str2double(get_field('kP'));
params.lineTension = str2double(get_field('line_tension'));
params.forceScale = str2double(get_field('force_scale')); % sigma_0, the force-scale!
params.connectivity = get_field('connectivity');
params.timeStep = str2double(get_field('time_step'));
params.dragCoeff = str2double(get_field('drag_coeff'));
params.jitterSize = l*str2double(get_field('jitter_std'));
% params.t1Threshold = 2;
params.t1Threshold = str2double(get_field('t1Threshold'));

% Grab the integration method
eval(['initialization.solver_method = ' ...
    get_field('solver_method') ';']);

% Set contractility gradient
% Activate "ventral fate"
eval(['x0 = ' get_field('x0') ';']);
eval(['xf = ' get_field('xf') ';']);
eval(['y0 = ' get_field('y0') ';']);
eval(['yf = ' get_field('yf') ';']);
contract.ventral.box = [ x0, y0, xf, yf];
contract.ventral.alt_tension = str2double(get_field('alt_tension'));

% Set the value of contractility in each cell
eval(['contract.contractility.spatial_model = ' ...
    get_field('contractModel') ';']);
eval(['contract.contractility.spatial_params = ' ...
    get_field('contractParams') ';']);
eval(['contract.contractility.temporal_model = ' ...
    get_field('tempModel') ';']);
eval(['contract.contractility.temporal_params = ' ...
    get_field('tempParams') ';']);

% Set the value of contractility in each cell
tis.setContractilityModel(contractions);

run_model(initialization,params,contract,OUT_DIR);
tisArr = assemble_model(OUT_DIR);

end