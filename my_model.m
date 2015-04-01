%%% energy model of ventral furrow formation
% close all; clear all; clear java;
% disp('Initializing...');

%% Set parameters
HEX_ANGLE = 'horizontal';
%HEX_ANGLE = 'vertical';
% HEX_ANGLE = 'diagonal';

HEX_NUM_X = 48;
HEX_NUM_Y = 16;

%Approx run times for different dimentions for 4 steps
% 48 by 16 - 70s initialize; 20s / step
% 36 by 12 - 12s initialize; 6 sec / step
% 24 by 8 - 10s initalize; 5 sec / step
% 12 by 4 - 5s initialize; 0.5 sec / step
l = 512 / max(HEX_NUM_X,HEX_NUM_Y) / 2;

TARGET_AREA_FRACTION_OF_INIT = 1;
AREA_ELASTICITY = 1e-4;
PERIM_ELASTICITY = 1e-2;
LINE_TENSION = 2;

CONNECTIVITY = 'purse string';
%  CONNECTIVITY = 'network';
% CONNECTIVITY = 'purse string and network';

STRATEGY = 'synchronous';
%   STRATEGY = 'random asynchronous';
% STRATEGY = 'neighbor asynchronous';

ELASTICITY = 'elastic';
% ELASTICITY = 'reorganization';

STEPS = 25; % number of constriction steps
TIME_STEP = 1;
VISCOSITY_COEFF = 10;

JITTERING_STD = l/5;

%% Initialize model

% create hexagons and the CellGraph object
tic
hexagons = create_hexagons(HEX_ANGLE,HEX_NUM_X, HEX_NUM_Y);
[centroid_list,regions] = get_cents(hexagons);
[vertex_list] = get_vertices(hexagons);
tis = Tissue(regions,vertex_list,centroid_list);

clear tissueArray; clf;
tissueArray(1:STEPS+1) = Tissue;
T = toc;
display(['Tissue initialized in ' num2str(T) ' sec'])
verts = tis.vert_coords;

A0 = mean([tis.getCells.area]);
P0 = mean([tis.getCells.perimeter]);

% Set parameters and initialize distance matrix
tic;
p.targetAreas = A0;
p.targetPerimeters = 0;
p.lineTension = LINE_TENSION;
p.areaElasticity = AREA_ELASTICITY;
p.perimElasticity = PERIM_ELASTICITY;
p.conn_opt = CONNECTIVITY;
p.viscosity = VISCOSITY_COEFF;
p.jitter_std = JITTERING_STD;
tis = tis.setParameters(p);
T = toc;
display(['Parameter and connection matrices initialized in ' num2str(T) ' sec'])

% Checks for parameter settings:
PERIM_ELASTICITY / AREA_ELASTICITY / A0 %Should be ~0.04
LINE_TENSION / AREA_ELASTICITY / A0^(3/2) / 2 %Should be ~0.12

%% JITTER + Show initial conditions
tic
tis_init = tis.jitterVertices(JITTERING_STD);
T = toc;
display(['Jitter added in ' num2str(T) ' sec'])
tis_init.draw(); title('Initial condition')

%% Euler scheme of model integration
% @todo: not numerically stable!

tis = tis_init; tissueArray(1) = tis; E = zeros(1,STEPS);

C = zeros( 1, tis.cells.length );
% C( [6] ) = 2e-4;
tis = tis.setContractility( C );

for i = 1:STEPS
    
    tic
    verts = tis.vert_coords;
    displacements = tis.get_force/p.viscosity;
    verts = verts + displacements;
%     keyboard
    tis = tis.evolve( verts );
    tissueArray( i + 1 ) = tis;
    E(i) = tis.get_energy;
    
    tis.draw('showVectors',displacements,'showActive'); drawnow;
    title(['Time step = ' num2str(i)])
    T = toc;
    display([num2str(i) '-th time step (' num2str(T) ' sec)'])
    
end

