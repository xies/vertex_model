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

AREA_ELASTICITY = 6e-5;
PERIM_ELASTICITY = 1e-4;
LINE_TENSION = 1;
FORCE_SCALE = 1; % sigma_0, the force-scale!

CONNECTIVITY = 'purse string';
%  CONNECTIVITY = 'network';
% CONNECTIVITY = 'purse string and network';

STEPS = 1000; % number of constriction steps
abs_tol = 1e-2; rel_tol = 1e-9;
TIME_STEP = 1e-8;
VISCOSITY_COEFF = 1e0;

JITTERING_STD = 1/10;

%% Initialize model

% create hexagons and the CellGraph object
tic
tis = Tissue(regions,vertex_list,centroid_list,CONNECTIVITY);

clear tissueArray p;
tissueArray(1:STEPS+1) = Tissue;
T = toc;
display(['Tissue initialized in ' num2str(T) ' sec'])
verts = tis.vert_coords;

A0 = mean([tis.getCells.area]);
P0 = mean([tis.getCells.perimeter]);
l = P0/6; % lattice length_scale
um_per_px = sqrt(40/A0); % pixel size

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
        't1Threshold',2 ...
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
        't1Threshold',2 ...
        };
end

tis.setParameters(param_config{:}); % set parameters
T = toc;
display(['Parameter and connection matrices initialized in ' num2str(T) ' sec'])

%% Set contractility gradient

tic

MODEL_FUN = @variable_cutoff;
CONTRACTILITY_MAGNITUDE = tis.parameters.areaElasticity*2;
CONT_STD = CONTRACTILITY_MAGNITUDE * 0.05;
CONTRACTILE_WIDTH = 40; % pxs
ALT_TENSION = 1;

% Activate "ventral fate"
box = [ 1/2-1/3 , 1/10 , 1/2+1/3 , 9/10 ];
cIDs = tis.getCellsWithinRegion(box);
tis.activateCell(cIDs,ALT_TENSION);
figure(1),tis.draw('showActive'); title('Ventral fated cells')

% Set the value of contractility in each cell
midline_x = tis.Xs/2; midline_y = tis.Ys/2;
contract_params = [CONTRACTILITY_MAGNITUDE , CONT_STD];
% contract_params = [CONTRACTILITY_MAGNITUDE midline_x, midline_y,...
%     CONTRACTILE_WIDTH CONT_STD];
tis.setContractilityModel(MODEL_FUN,contract_params);
C = tis.getContractility;
T = toc;
display(['Contractility set in ' num2str(T) ' sec'])
figure(2),tis.draw('showContractile'); title('Contractility')

%% JITTER + Show initial conditions

% Add some jitter
tic
tis.jitterVertices(JITTERING_STD);
num_cells = tis.cells.length;

T = toc;
display(['Jitter added and contractility set in ' num2str(T) ' sec'])
tis.draw('showContractile'); title('Initial condition')

%% Runge-Kutta 2/3

opt = odeset('OutputFcn',@odeprint);
OUT_DIR = '~/Desktop/tmp';

%% Euler scheme of model integration

tissueArray(1) = tis;
E = zeros(1,STEPS); Econt = zeros(1,STEPS);
contractility = zeros(num_cells,STEPS);

for i = 1:100
    
    tic
    
    verts_old = tis.vert_coords;
    displacements = tis.get_force ...
        /tis.parameters.viscosity * tis.parameters.lengthScale ...
        * tis.parameters.stepSize / um_per_px;
    verts = verts_old + displacements;
    
    tis = tis.evolve( verts );
    tissueArray( i + 1 ) = tis;
    E(i) = tis.get_energy;
    
    if (abs(E - E_prev) < rel_tol * E ...
            || abs(E - E_prev) < abs_tol), break; end
    while i>1 && E(i) - E(i-1) > 0
        % If unstable, try again with smaller stepSize (by 1/5) each time
%         display(['New step size: ' num2str(tis.parameters.stepSize / 5)])
%         tis_old = tissueArray(i);
%         tis_old.parameters.stepSize = tis.parameters.stepSize / 5;
%         verts = tis_old.vert_coords;
%         displacements = tis_old.get_force ...
%             / tis.parameters.viscosity * tis.parameters.lengthScale ...
%             * tis.parameters.stepSize / um_per_px;
%         verts = verts + displacements;
%         
%         tis = tis_old.evolve( verts );
%         tis.parameters.stepSize = tis.parameters.stepSize/5;
%         tissueArray( i + 1 ) = tis;
%         E(i) = tis.get_energy;
%         keyboard
        error('Unstable!')
    end
    
%     if mod(i,10) == 0
%         % Add random noise
%         C(cIDs) = max(C(cIDs) + CONTRACTILITY_MAGNITUDE/100 * randn(1,numel(cIDs)),0);
% %         C(cIDs) = min(C(cIDs),5e-4);
%         tis = tis.setContractility(C);
%     end
    
%     figure(1)
%     tis.draw('showVectors',displacements, ...
%         'showContractile');
%     title(['Time step = ' num2str(i)]);
%     figure(2)
%     hist(displacements(:),30);
%     drawnow;
    
    T = toc;
    display([num2str(i) '-th time step (' num2str(T) ' sec)'])
    
end

tissueArray(i+1:end) = [];
% save(['~/Dropbox (MIT)/model_contractility_to_Ka_1' ],'tissueArray');
