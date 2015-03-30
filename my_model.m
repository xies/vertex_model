%%% energy model of ventral furrow formation
% close all; clear all; clear java;
% disp('Initializing...');

%% set parameters
HEX_ANGLE = 'horizontal';
%HEX_ANGLE = 'vertical';
% HEX_ANGLE = 'diagonal';

HEX_NUM_X = 4;
HEX_NUM_Y = 5;

%Approx run times for different dimentions for 4 steps
%48 by 16 - 26 min
%36 by 12 - 5.2 min
%24 by 8 - 23 sec
%12 by 4 - 5 sec

TARGET_AREA_FRACTION_OF_INIT = 1;
AREA_ELASTICITY = 1e-4;
PERIM_ELASTICITY = 4e-2;
LINE_TENSION = 20;

CONNECTIVITY = 'purse string';
%  CONNECTIVITY = 'network';
% CONNECTIVITY = 'purse string and network';

STRATEGY = 'synchronous';
%   STRATEGY = 'random asynchronous';
% STRATEGY = 'neighbor asynchronous';

ELASTICITY = 'elastic';
% ELASTICITY = 'reorganization';

STEPS = 50; % number of constriction steps
TIME_STEP = 1;
VISCOSITY_COEFF = 5;

JITTERING_STD = 10;

%% Create hexagons and the CellGraph object

hexagons = create_hexagons(HEX_ANGLE,HEX_NUM_X, HEX_NUM_Y);
[centroid_list,regions] = get_cents(hexagons);
[vertex_list] = get_vertices(hexagons);

tis = Tissue(regions,vertex_list,centroid_list);
verts = tis.vert_coords;

A0 = mean([tis.getCells.area]);
P0 = mean([tis.getCells.perimeter]);

% Set parameters and initialize distance matrix
p.targetAreas = A0;
p.targetPerimeters = 0;
p.lineTension = LINE_TENSION;
p.areaElasticity = AREA_ELASTICITY;
p.perimElasticity = PERIM_ELASTICITY;
p.conn_opt = CONNECTIVITY;
p.viscosity = VISCOSITY_COEFF;
tis = tis.setParameters(p);

% Checks for parameter settings:
PERIM_ELASTICITY / AREA_ELASTICITY / A0 %Should be ~0.04
LINE_TENSION / AREA_ELASTICITY / A0^(3/2) / 2 %Should be ~0.12

% JITTER + Show initial conditions
clear tissueArray;
tissueArray(1:STEPS+1) = Tissue;
tis_init = tis.jitterVertices(JITTERING_STD);
tis_init.draw(); title('Initial condition')

%% Euler scheme of model integration
% @todo: not numerically stable!

tis = tis_init; tissueArray(1) = tis; E = zeros(1,STEPS);

C = zeros( 1, tis.cells.length );
% C( [6] ) = 2e-4;
tis = tis.setContractility( C );

for i = 1:STEPS
    
    verts = tis.vert_coords;
    displacements = tis.get_force/p.viscosity;
    verts = verts + displacements;
%     keyboard
    tis = tis.evolve( verts );
    tissueArray( i + 1 ) = tis;
    E(i) = tis.get_energy
    
    tis.draw('showVectors',displacements,'showActive'); drawnow;
    title(['Time step = ' num2str(i)])
    i
    
end

