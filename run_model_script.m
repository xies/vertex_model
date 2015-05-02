%%% Vertex model of ventral furrow formation

%Approx run times for different dimentions for 4 steps
% 48 by 16 - 70s initialize; 20s / step
% 36 by 12 - 12s initialize; 6 sec / step
% 24 by 8 - 10s initalize; 5 sec / step
% 12 by 4 - 5s initialize; 0.5 sec / step

%% Set parameters

initialization.initialize = @create_hexagons;
numX = 8; numY = 8;

initialization.model_params = {'diagonal',numX,numY};
initialization.steps = 5000; % number of constriction steps
l = 512 / max(numX,numY) / 2;

params.dimensionless = 0;
params.areaElasticity = 2.5e-12;
params.perimElasticity = 0;
params.lineTension = 0.2e-12;
params.forceScale = 1; % sigma_0, the force-scale!
params.connectivity = 'purse string';
params.timeStep = 0.01;
params.dragCoeff = 1e3;
params.jitterSize = l/5;

%% Set contractility gradient

tic
% Activate "ventral fate"
contract.ventral.box = [ Inf, Inf ...
    Inf, Inf];
contract.ventral.alt_tension = 1;

% Set the value of contractility in each cell
contract.contractility.model = @uniform_cutoff;
contractMag = 500; % times areaElasticity
contractVar = contractMag * params.areaElasticity;
% gradientWidth = 40; % pxs
% midline_x = tis_init.Xs/2; midline_y = tis_init.Ys/2;
contract.contractility.params = [contractMag, contractVar];

tisArr = run_model(initialization,params,contract);

% save()