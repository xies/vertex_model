%%% Vertex model of ventral furrow formation

%% Set parameters
clc
% HEX_ANGLE = 'horizontal';
% HEX_ANGLE = 'vertical';
HEX_ANGLE = 'diagonal';

HEX_NUM_X = 8;
HEX_NUM_Y = 8;
hexagons = create_hexagons(HEX_ANGLE,HEX_NUM_X, HEX_NUM_Y);
[centroid_list,regions] = get_cents(hexagons);
[vertex_list] = get_vertices(hexagons);

%Approx run times for different dimentions for 4 steps
% 48 by 16 - 70s initialize; 20s / step
% 36 by 12 - 12s initialize; 6 sec / step
% 24 by 8 - 10s initalize; 5 sec / step
% 12 by 4 - 5s initialize; 0.5 sec / step
l = max(size(regions)) / max(HEX_NUM_X,HEX_NUM_Y) / 2;

DIMENSIONLESS = 0;

AREA_ELASTICITY = 2.5e-12;
PERIM_ELASTICITY = 0;
LINE_TENSION = 0.2e-12;
FORCE_SCALE = 1; % sigma_0, the force-scale!

CONNECTIVITY = 'purse string';
%  CONNECTIVITY = 'network';
% CONNECTIVITY = 'purse string and network';

STEPS = 10000; % number of constriction steps
TIME_STEP = 0.01;
VISCOSITY_COEFF = 1e3;

JITTERING_STD = l/5;

%% Initialize model

% create hexagons and the CellGraph object
tic
tis_init = Tissue(regions,vertex_list,centroid_list,CONNECTIVITY);

clear tissueArray p;
tissueArray(1:STEPS+1) = Tissue;
T = toc;
display(['Tissue initialized in ' num2str(T) ' sec'])
verts = tis_init.vert_coords;

A0 = mean([tis_init.getCells.area]);
P0 = mean([tis_init.getCells.perimeter]);
l = P0/6; % lattice length_scale
um_per_px = sqrt(0.02/A0); % pixel size

if DIMENSIONLESS
    param_config = {...
        'dimensonless', true, ...
        'lengthScale', l, ...
        'forceScale', FORCE_SCALE, ...
        'lineTension', LINE_TENSION, ...
        'lineAnisotropy', 1, ...
        'areaElasticity', AREA_ELASTICITY*l^3/FORCE_SCALE, ...
        'perimElasticity', PERIM_ELASTICITY*l/FORCE_SCALE, ...
        'targetArea',sqrt(3)*3/2, ...
        'viscosity', VISCOSITY_COEFF * FORCE_SCALE / l, ...
        'connectivity', CONNECTIVITY, ...
        'stepSize', TIME_STEP, ...
        'jitterSize', JITTERING_STD, ...
        'um_per_px', um_per_px, ...
        'dt_per_frame', 10, ...
        't1Threshold',3 ...
        };
else
    param_config = {...
        'dimensonless', false, ...
        'lengthScale', 1, ...
        'forceScale', FORCE_SCALE, ...
        'lineTension', LINE_TENSION, ...
        'lineAnisotropy', 1, ...
        'areaElasticity', AREA_ELASTICITY, ...
        'perimElasticity', PERIM_ELASTICITY, ...
        'viscosity', VISCOSITY_COEFF, ...
        'connectivity', CONNECTIVITY, ...
        'stepSize', TIME_STEP, ...
        'jitterSize', JITTERING_STD, ...
        'um_per_px', um_per_px, ...
        'dt_per_frame', 10 ...
        't1Threshold',0 ...
        };
end

tis_init = tis_init.setParameters(param_config{:}); % set parameters
T = toc;
display(['Parameter and connection matrices initialized in ' num2str(T) ' sec'])

%% Set contractility gradient

tic

MODEL_FUN = @uniform_cutoff;
CONTRACTILITY_MAGNITUDE = tis_init.parameters.areaElasticity*500;
CONT_STD = CONTRACTILITY_MAGNITUDE * .2;
CONTRACTILE_WIDTH = 40; % pxs
ALT_TENSION = tis_init.parameters.areaElasticity*1;

% Activate "ventral fate"
box = [ tis_init.Xs * (1/2-1/12), tis_init.Ys * (1/2-1/12) ...
    tis_init.Xs * (1/2+1/12), tis_init.Ys * (1/2+1/12)];
cIDs = tis_init.getCellsWithinRegion(box);
tis_init = tis_init.activateCell(cIDs,ALT_TENSION);
figure(1),tis_init.draw('showActive'); title('Ventral fated cells')

% Set the value of contractility in each cell
midline_x = tis_init.Xs/2; midline_y = tis_init.Ys/2;
contract_params = [CONTRACTILITY_MAGNITUDE, CONT_STD];
% contract_params = [CONTRACTILITY_MAGNITUDE midline_y,...
%     CONTRACTILE_WIDTH];
tis_init = tis_init.setContractilityModel(MODEL_FUN,contract_params);
C = tis_init.getContractility;
T = toc;
display(['Contractility set in ' num2str(T) ' sec'])
figure(2),tis_init.draw('showContractile'); title('Contractility')

%% JITTER + Show initial conditions

% Add some jitter
tic
tis_init = tis_init.jitterVertices(JITTERING_STD);
num_cells = tis_init.cells.length;

T = toc;
display(['Jitter added and contractility set in ' num2str(T) ' sec'])
tis_init.draw('showContractile'); title('Initial condition')

%% Euler scheme of model integration

tis = tis_init; tissueArray(1) = tis;
E = zeros(1,STEPS); Econt = zeros(1,STEPS);
contractility = zeros(num_cells,STEPS);

for i = 1:STEPS
    
    tic
    
    verts = tis.vert_coords;
    displacements = tis.get_force ...
        /tis.parameters.viscosity * tis.parameters.lengthScale ...
        * tis.parameters.stepSize;
    verts = verts + displacements;
    
    tis = tis.evolve( verts );
    tissueArray( i + 1 ) = tis;
    E(i) = tis.get_energy;
    
    if i>1 && abs(E(i) - E(i-1)) < 1e-5, break; end
    if i>1 && E(i) - E(i-1) > 0, error('Unstable regime!'); end
    
%     if mod(i,10) == 0
%         % Add random noise
%         C(cIDs) = max(C(cIDs) + CONTRACTILITY_MAGNITUDE/100 * randn(1,numel(cIDs)),0);
% %         C(cIDs) = min(C(cIDs),5e-4);
%         tis = tis.setContractility(C);
%     end
    
%     tis.draw('showVectors',displacements, ...
%         'showContractile');
%     title(['Time step = ' num2str(i)]);
%     figure(2)
%     hist(displacements(:),30);
%     figure(1)
%     drawnow;
    
    T = toc;
    display([num2str(i) '-th time step (' num2str(T) ' sec)'])
    
end

tissueArray(i+1:end) = [];
% save(['~/Dropbox (MIT)/model_contractility_to_Ka_1' ],'tissueArray');