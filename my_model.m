%%% energy model of ventral furrow formation
% close all; clear all; clear java;
% disp('Initializing...');

%% set parameters
HEX_ANGLE = 'horizontal';
%HEX_ANGLE = 'vertical';
% HEX_ANGLE = 'diagonal';

HEX_NUM_X = 16;
HEX_NUM_Y = 8;

%Approx run times for different dimentions for 4 steps
%48 by 16 - 26 min
%36 by 12 - 5.2 min
%24 by 8 - 23 sec
%12 by 4 - 5 sec

TARGET_AREA_FRACTION_OF_INIT = 1;
AREA_ELASTICITY = 1;

CONNECTIVITY = 'purse string';
%  CONNECTIVITY = 'network';
% CONNECTIVITY = 'purse string and network';

STRATEGY = 'synchronous';
%   STRATEGY = 'random asynchronous';
% STRATEGY = 'neighbor asynchronous';

ELASTICITY = 'elastic';
% ELASTICITY = 'reorganization';

STEPS = 50; % number of constriction steps

PASSIVE_LAYER_THICKNESS = 1;
JITTERING_STD = 1;

INITIAL_EQM_LENGTH_FACTOR = 1/2;

FRACTION_OF_ACTIVE_CELLS_TO_CONSTRICT = 0.3;

%% create hexagons and the CellGraph object

hexagons = create_hexagons(HEX_ANGLE, HEX_NUM_X, HEX_NUM_Y);
[centroid_list,regions] = get_cents(hexagons);
[vertex_list] = get_vertices(hexagons);

tis = Tissue(regions,vertex_list,centroid_list);
verts = tis.vert_coords;

%% Set parameters

A0 = mean([tis.getCells.area]);
P0 = mean([tis.getCells.perimeter]);

p.targetAreas = A0;
p.targetPerimeter = P0;
p.areaElasticity = AREA_ELASTICITY;
p.conn_opt = CONNECTIVITY;

tis = tis.setParameters(p);

clear tissueArray;
tissueArray(1:STEPS+1) = Tissue;
tis_init = tis;

%% Euler scheme of model integration
% @todo: not numerically stable!

tis = tis_init; tissueArray(1) = tis;

for i = 1:STEPS
    
    tis = tis.jitterVertices( JITTERING_STD );
    verts = tis.vert_coords;
    displacements = tis.get_velocities / 1e5;
    verts = verts + displacements;
    
    tis = tis.evolve( verts );
    tissueArray( i + 1 ) = tis;
    i
    
end

