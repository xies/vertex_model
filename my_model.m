%%% energy model of ventral furrow formation
% close all; clear all; clear java;
% disp('Initializing...');

%% set parameters
HEX_ANGLE = 'horizontal';
%HEX_ANGLE = 'vertical';
% HEX_ANGLE = 'diagonal';

HEX_NUM_X = 36;
HEX_NUM_Y = 12;

%Approx run times for different dimentions for 4 steps
%48 by 16 - 26 min
%36 by 12 - 5.2 min
%24 by 8 - 23 sec
%12 by 4 - 5 sec

SPRING_CONSTANT_INITIAL = 1;

CONNECTIVITY = 'purse string';
%  CONNECTIVITY = 'network';
% CONNECTIVITY = 'purse string and network';

STRATEGY = 'synchronous';
%   STRATEGY = 'random asynchronous';
% STRATEGY = 'neighbor asynchronous';

ELASTICITY = 'elastic';
% ELASTICITY = 'reorganization';

STEPS = 4; % number of constriction steps
MEAN_NUMBER_OF_CONSTRICTIONS = 2;

PASSIVE_LAYER_THICKNESS = 1;
JITTERING_STD = 1;
FRACTIONAL_NEW_EQUILIBRIUM_LENGTH = .5;
FRACTIONAL_NEW_SPRING_CONSTANT = 1.1;

INITIAL_EQM_LENGTH_FACTOR = 1/2;

FRACTION_OF_ACTIVE_CELLS_TO_CONSTRICT = 0.5;

%% create hexagons and the CellGraph object
hexagons = create_hexagons(HEX_ANGLE, HEX_NUM_X, HEX_NUM_Y);
[centroid_list,regions] = get_cents(hexagons);
[vertex_list] = get_vertices(hexagons);

tis = Tissue(regions,vertex_list,centroid_list);
verts = tis.vert_coords;

%% build the parameter cell array
% fixed cells have at most 2 neighbors

fixed_cells_logical = tis.numCellTouchingVertices <= 2;

p.fixed_cells = find(fixed_cells_logical);
p.not_fixed_cells = find(~fixed_cells_logical);
p.numV = size(verts, 1);

p.preferred_distances = INITIAL_EQM_LENGTH_FACTOR * squareform( pdist(verts) );
%This is where you can change the spring constant
p.spring_constants = SPRING_CONSTANT_INITIAL*ones(size(p.preferred_distances));

clear tissueArray;
tissueArray(1:STEPS+1) = Tissue;
P = cell(STEPS, 1);
Force = cell(STEPS, 3); % 3 rows -> bigvert, biggrad, stress

r.average_anisotropy = zeros(STEPS, 1);
r.average_xy_anisotropy = zeros(STEPS, 1);
r.average_constriction = zeros(STEPS, 1);
r.average_constriction_std = zeros(STEPS, 1);
r.rms_stress = zeros(STEPS, 1);
r.max_stress = zeros(STEPS, 1);
r.stress_anisotropy = zeros(STEPS, 1);

%% initial set-up
% Make outer layer passive
% cg.setActiveCellsAuto(PASSIVE_LAYER_THICKNESS);
tis = tis.connectVertices(CONNECTIVITY);
p.connectivity = tis.connectivity;

% % so far connectivity doesn't change each round
% % but with the stress reponse fibers it will
% switch CONNECTIVITY
%     case 'purse string'
%         p.connectivity = cg.connectivityMatrixVertex;
%     case 'network'
%         p.connectivity = cg.connectivityMatrixVertexCrossConnections;
%     case 'purse string and network'
%         p.connectivity = cg.connectivityMatrixVertex + ...
%             cg.connectivityMatrixVertexCrossConnections;
%     otherwise
%         disp('No valid CONNECTIVITY selected');
% end

%% Minimize energy
tissueArray(1) = tis;
disp('Beginning constriction...');

for j = 1:STEPS
    
    display(['Time step = ' num2str(j)])
    
    % Add jitter
    verts = verts + JITTERING_STD*randn(size(verts));
    
    tis = tis.deactivateCell;
    tis = tis.activateCell( 'random', FRACTION_OF_ACTIVE_CELLS_TO_CONSTRICT );
    tis = tis.deactivateBorder;
    
    if isempty(tis.getActiveCells)
        actconn = zeros(size(p.connectivity));
    else
        actconn = tis.connectVertices( 'purse string', tis.getActiveCells ).connectivity;
    end
    actconn = logical(actconn);
    
    % choose cells to constrict out of the active cells
    p.preferred_distances(actconn) = ...
        p.preferred_distances(actconn) * FRACTIONAL_NEW_EQUILIBRIUM_LENGTH;
    p.spring_constants(actconn) = ...
        p.spring_constants(actconn) * FRACTIONAL_NEW_SPRING_CONSTANT;
    
    p.initial_verts = verts;
    
    % optimization routine
    % use the anonymous function @(x) calc_energy(x,p) to pass parameters...
    options = optimset('GradObj', 'on', 'DerivativeCheck', 'off', ...
        'LargeScale', 'off', 'Hessian', 'off', 'display','off');
    
    % remove the fixed cells from "verts"
    verts(p.fixed_cells, :) = [];
    % change the shape of verts to be a column vector (note: y before x)
    verts = [verts(:,1); verts(:,2)];
    
    [output_verts, E, exitflag, output, grad_true] ...
        = fminunc(@(x) calculate_energy(x, p), verts, options);
    output_verts = reshape(output_verts, length(output_verts)/2, 2);  % change it back
    % put the original vertices back in the fixed cells
    verts = p.initial_verts;
    verts(p.not_fixed_cells, :) = output_verts;
    
    % update tissue with new vertex coordinates
    tis = tis.evolve(verts);
    tissueArray(j+1) = tis;
    
    % store the new CellGraph
    %     cgArray(j) = CellGraph(cg);
    % reset the active cells (if desired???)
    %     cg.setActiveCellsAuto(PASSIVE_LAYER_THICKNESS)
    
end
disp('Constriction finished.');

%% display results

average_anisotropy = r.average_anisotropy
average_xy_anisotropy = r.average_xy_anisotropy
average_constriction = r.average_constriction
average_constriction_std = r.average_constriction_std
rms_stress = r.rms_stress
max_stress = r.max_stress
stress_anisotropy = r.stress_anisotropy
stress_anis = r.stress_anisotropy

% plotting
drawing(:,:,3) = initial_cg.draw;
drawing(:,:,1) = initial_cg.drawVertices;
figure('Position', [100 460 600 300]);
imshow(drawing);
title('Beginning graph')

%%

for j = STEPS:STEPS
    
    cg = cgArray(j);
    verts = cg.vertexCoords;
    p = P{j};
    bigvert = Force{j, 1};
    biggrad = Force{j, 2};
    stress  = Force{j, 3};
    
    drawing(:,:,3) = cg.draw;
    drawing(:,:,1) = cg.drawVertices;
    figure('Position', [100 50 600 400]);
    imshow(drawing);
end

for j = STEPS:STEPS
    figure('Position', [800 50 600 400]);
    imshow(drawing);
    hold on;
    
    % plot active cells in yellow
    cents = cg.centroidCoordsActiveCells;
    plot(cents(:,2), cents(:,1), '*y');
    
    quiver(bigvert(:, 2), bigvert(:,1), ...
        biggrad(:, 2), biggrad(:, 1), 0.5, 'g');
    hold off
end
for j = STEPS:STEPS
    
    stress_img = zeros(size(cg.draw));
    for k = 1:size(verts, 1)
        y = round(verts(k, 1)); x = round(verts(k, 2));
        s = 2;
        stress_img(y-s:y+s, x-s:x+s) = stress(k);
    end
    figure('Position', [800 460 600 300]); hold on;
    imshow(zeros(size(stress_img))); %to get the axes right
    imagesc(stress_img); axis equal; axis off;
    colorbar; colormap 'hot'; zoom(2); zoom(1/2);
    draw_connectivity_matrix_fast(cg.connectivityMatrixVertex, ...
        verts, [0 0 1]);
    hold off;
    
end

%% save p
save('Ps_for_all_steps.mat','P','cg','initial_cg')
