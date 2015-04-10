%%% Vertex model of ventral furrow formation

%% Set parameters
clc
% HEX_ANGLE = 'horizontal';
% HEX_ANGLE = 'vertical';
HEX_ANGLE = 'diagonal';

HEX_NUM_X = 10;
HEX_NUM_Y = 10;
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

AREA_ELASTICITY = 2.3e-1;
PERIM_ELASTICITY = 2e-2;
FORCE_SCALE = 1; % sigma_0, the force-scale!

CONNECTIVITY = 'purse string';
%  CONNECTIVITY = 'network';
% CONNECTIVITY = 'purse string and network';

STRATEGY = 'synchronous';
%   STRATEGY = 'random asynchronous';
% STRATEGY = 'neighbor asynchronous';

ELASTICITY = 'elastic';
% ELASTICITY = 'reorganization';

STEPS = 1000; % number of constriction steps
TIME_STEP = 0.01;
VISCOSITY_COEFF = 1;

JITTERING_STD = l/10;

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
um_per_px = sqrt(40/A0); % pixel size

if DIMENSIONLESS
    param_config = {...
        'dimensonless', true, ...
        'lengthScale', l, ...
        'forceScale', FORCE_SCALE, ...
        'lineTension', FORCE_SCALE, ...
        'lineAnisotropy', 2, ...
        'areaElasticity', AREA_ELASTICITY*l^3/FORCE_SCALE, ...
        'perimElasticity', PERIM_ELASTICITY*l/FORCE_SCALE, ...
        'targetArea',sqrt(3)*3/2, ...
        'viscosity', VISCOSITY_COEFF * FORCE_SCALE / l, ...
        'connectivity', CONNECTIVITY, ...
        'stepSize', TIME_STEP, ...
        'jitterSize', JITTERING_STD, ...
        'um_per_px', um_per_px, ...
        'dt_per_frame', 10, ...
        't1Threshold',1 ...
        };
else
    param_config = {...
        'dimensonless', false, ...
        'lengthScale', 1, ...
        'forceScale', FORCE_SCALE, ...
        'lineTension', FORCE_SCALE, ...
        'lineAnisotropy', 2, ...
        'areaElasticity', AREA_ELASTICITY, ...
        'perimElasticity', PERIM_ELASTICITY, ...
        'viscosity', VISCOSITY_COEFF, ...
        'connectivity', CONNECTIVITY, ...
        'stepSize', TIME_STEP, ...
        'jitterSize', JITTERING_STD, ...
        'um_per_px', um_per_px, ...
        'dt_per_frame', 10 ...
        't1Threshold',1 ...
        };
end

tis_init = tis_init.setParameters(param_config{:}); % set parameters
T = toc;
display(['Parameter and connection matrices initialized in ' num2str(T) ' sec'])

%% Set contractility gradient
tic

% Activate "ventral fate"
box = [ tis_init.Xs* 1/4, tis_init.Ys * 1/10 ...
    tis_init.Xs * 3/4, tis_init.Ys * 9/10];
cIDs = tis_init.getCellsWithinRegion(box);
tis_init = tis_init.activateCell(cIDs);
figure(1),tis_init.draw('showActive'); title('Ventral fated cells')

tis_init = tis_init.setContractilityModel()
T = toc;
display(['Contractility set in ' num2str(T) ' sec'])
figure(1),tis_init.draw('showActive'); title('Contractility')

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
    
    cy = tis_init.getCentroidY;
%     C( ind ) = max( C(ind) + AREA_ELASTICITY* (5*randn(1,numel(ind))+5), 0);
    tis = tis.activateCell( find(C > 0) );
    tis = tis.setContractility( C );
%     if mod(i,10) == 0,
%         bias = AREA_ELASTICITY*( 0.05*randn(1,numel(ind)) );
%     end
    
    verts = tis.vert_coords;
    displacements = tis.get_force/p.viscosity * p.length_scale * p.step_size;
    verts = verts + displacements;
    
    if max(displacements) < .05,
        break;
    end
    
    tis = tis.evolve( verts );
    tissueArray( i + 1 ) = tis;
    E(i) = tis.get_energy;
    
    tis.draw('showVectors', p.dt* displacements,'showContractile');
    title(['Time step = ' num2str(i)]);
    figure(2)
    hist(displacements(:),30);
    figure(1)
    drawnow;
    
    T = toc;
    display([num2str(i) '-th time step (' num2str(T) ' sec)'])
    
end

%% Use MATLAB's RK45 ODE solver (ODE45)

% ode_solver = @ode45;
initial_conditions = tis_init.vert_coords;
tis = tis_init;
% 
opt = odeset('InitialStep',0.01);
[T,Y] = ode_solver( @(t,y) evolve_tissue(t,y,tis), ...
    [0 60],initial_conditions,opt );

