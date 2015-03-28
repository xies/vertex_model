%%% energy model of ventral furrow formation
% close all; clear all; clear java;
% disp('Initializing...');

%% set parameters
HEX_ANGLE = 'horizontal';
%HEX_ANGLE = 'vertical';
% HEX_ANGLE = 'diagonal';

HEX_NUM_X = 8;
HEX_NUM_Y = 6;

%Approx run times for different dimentions for 4 steps
%48 by 16 - 26 min
%36 by 12 - 5.2 min
%24 by 8 - 23 sec
%12 by 4 - 5 sec

TARGET_AREA_FRACTION_OF_INIT = 1;
AREA_ELASTICITY = 1e-4;
PERIM_ELASTICITY = .01;
LINE_TENSION = 10;

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

JITTERING_STD = 2;

%% Create hexagons and the CellGraph object

hexagons = create_hexagons(HEX_ANGLE,HEX_NUM_X, HEX_NUM_Y);
[centroid_list,regions] = get_cents(hexagons);
[vertex_list] = get_vertices(hexagons);

tis = Tissue(regions,vertex_list,centroid_list);
verts = tis.vert_coords;

A0 = mean([tis.getCells.area]);
P0 = mean([tis.getCells.perimeter]);

p.targetAreas = A0;
p.targetPerimeters = P0;
p.lineTension = LINE_TENSION;
p.areaElasticity = AREA_ELASTICITY;
p.perimElasticity = PERIM_ELASTICITY;
p.conn_opt = CONNECTIVITY;

tis = tis.setParameters(p);

clear tissueArray;
tissueArray(1:STEPS+1) = Tissue;
tis_init = tis;

%% Euler scheme of model integration
% @todo: not numerically stable!

tis = tis_init; tissueArray(1) = tis; E = zeros(1,STEPS);
tis = tis.jitterVertices(JITTERING_STD);

C = zeros( 1, tis.cells.length );
C( [8] ) = 0;
tis = tis.setContractility( C );

for i = 1:STEPS
    
    verts = tis.vert_coords;
    displacements = tis.get_velocities/1;
    verts = verts + displacements;
    
    
    I = tis.draw(); imagesc(I); axis square, hold on;
    scatter(tis.vert_coords(:,2),tis.vert_coords(:,1),100,'w')
    quiver(tis.vert_coords(:,2),tis.vert_coords(:,1),displacements(:,2),displacements(:,1),0,'w-');
    drawnow;
    
    tis = tis.evolve( verts );
    tissueArray( i + 1 ) = tis;
    E(i) = tis.get_energy
    
    i
%     keyboard
    
end

