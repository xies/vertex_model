%%% Vertex model of ventral furrow formation

%% Set parameters
% HEX_ANGLE = 'horizontal';
HEX_ANGLE = 'vertical';
% HEX_ANGLE = 'diagonal';

HEX_NUM_X = 25;
HEX_NUM_Y = 15;
hexagons = create_hexagons(HEX_ANGLE,HEX_NUM_X, HEX_NUM_Y);
[centroid_list,regions] = get_cents(hexagons);
[vertex_list] = get_vertices(hexagons);

%Approx run times for different dimentions for 4 steps
% 48 by 16 - 70s initialize; 20s / step
% 36 by 12 - 12s initialize; 6 sec / step
% 24 by 8 - 10s initalize; 5 sec / step
% 12 by 4 - 5s initialize; 0.5 sec / step
l = max(size(regions)) / max(HEX_NUM_X,HEX_NUM_Y) / 2;

TARGET_AREA_FRACTION_OF_INIT = 1;
AREA_ELASTICITY = 5e-4;
PERIM_ELASTICITY = 1e-2;
LINE_TENSION = 1; % sigma_0, the force-scale!

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
VISCOSITY_COEFF = 1e3;

JITTERING_STD = l/10;

%% Initialize model

% create hexagons and the CellGraph object
tic
tis_init = Tissue(regions,vertex_list,centroid_list);

clear tissueArray p;
tissueArray(1:STEPS+1) = Tissue;
T = toc;
display(['Tissue initialized in ' num2str(T) ' sec'])
verts = tis_init.vert_coords;

A0 = mean([tis_init.getCells.area]);
P0 = mean([tis_init.getCells.perimeter]);
l = P0/6; % lattice length_scale
p.um_per_px = sqrt(40/A0); % pixel size

% Set parameters and initialize distance matrix
% Use only dimension-less reduced parameters!
tic;
p.step_size = TIME_STEP; % Numerical method step size
p.dt = TIME_STEP*10; % Estimate time (in seconds)
p.length_scale = l;
p.targetAreas = sqrt(3)*3/2;
p.targetPerimeters = 0;
p.lineTension = 1;
p.areaElasticity = AREA_ELASTICITY*l^3/LINE_TENSION;
p.perimElasticity = PERIM_ELASTICITY*l/LINE_TENSION;
p.conn_opt = CONNECTIVITY;
p.viscosity = VISCOSITY_COEFF * LINE_TENSION / l;
p.jitter_std = JITTERING_STD; % px, not dimensionless
p.t1_threshold = 3; %px, not dimensionless!
tis_init = tis_init.setParameters(p); % set parameters
T = toc;
display(['Parameter and connection matrices initialized in ' num2str(T) ' sec'])

% Checks for parameter settings:
% PERIM_ELASTICITY / AREA_ELASTICITY / (A0 * um_per_px^2 ) %Should be ~0.04
% LINE_TENSION / AREA_ELASTICITY / (A0 * um_per_px^2 )^(3/2) / 2 %Should be ~0.12

%% JITTER + Show initial conditions

% Add some jitter
tic
tis_init = tis_init.jitterVertices(JITTERING_STD);
num_cells = tis_init.cells.length;

% Set contractility gradient
C = zeros( 1, num_cells );
ind = tis_init.getCellsWithinRegion( ...
    [tis_init.Xs*1/4, tis_init.Ys*1/10,...
    tis_init.Xs*3/4, tis_init.Ys*9/10] );
cy = tis_init.getCentroidY; midline = tis_init.Xs/2;
C( ind ) = 5*exp(-(cy(ind)-midline).^2 / 50^2);
tis_init = tis_init.activateCell( find(C > 0) );
tis_init = tis_init.setContractility( C );
tis_init = tis_init.deactivateBorder;

bias = 1e-2*randn(1,numel(ind));

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
    C( ind ) = max( C(ind) + 5e-2 * randn(1,numel(ind)) + bias, 0);
    tis = tis.activateCell( find(C > 0) );
    tis = tis.setContractility( C );
    
    verts = tis.vert_coords;
    displacements = tis.get_force/p.viscosity * p.length_scale * p.step_size;
    verts = verts + displacements;
    
    if max(displacements) < .01,
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

