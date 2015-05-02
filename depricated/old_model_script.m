%%% energy model of ventral furrow formation
close all; clear all; clear java; clear java; clear java; clear java; clc; clear java;
disp('Initializing...');

%% set parameters
 %HEX_ANGLE = 'horizontal';
 %HEX_ANGLE = 'vertical';
HEX_ANGLE = 'diagonal';


HEX_NUM_X = 12;
HEX_NUM_Y = 4;

%Approx times for different dimentions
%48 by 16 - 26 min
%36 by 12 - 5.2 min
%24 by 8 - 23 sec
%12 by 4 - 5 sec

SPRING_CONSTANT_INITIAL = 1;

CONNECTIVITY = 'purse string';
%  CONNECTIVITY = 'network';
% CONNECTIVITY = 'purse string and network';

% STRATEGY = 'synchronous';
  STRATEGY = 'random asynchronous';
% STRATEGY = 'neighbor asynchronous';

%ELASTICITY = 'elastic';
ELASTICITY = 'reorganization';

STEPS = 4;  % number of constriction steps
MEAN_NUMBER_OF_CONSTRICTIONS = 2;

PASSIVE_LAYER_THICKNESS = 0;
JITTERING_STD = 0;
FRACTIONAL_NEW_EQUILIBRIUM_LENGTH = .5;
FRACTIONAL_NEW_SPRING_CONSTANT = 2;

INITIAL_EQM_LENGTH_FACTOR = 1/2;

FRACTION_OF_ACTIVE_CELLS_TO_CONSTRICT = MEAN_NUMBER_OF_CONSTRICTIONS / STEPS;

%% create hexagons and the CellGraph object
hexagons = create_hexagons(HEX_ANGLE, HEX_NUM_X, HEX_NUM_Y);
[centroid_list regions] = get_cents(hexagons);
[vertex_list] = get_vertices(hexagons);
% below: need to make integers so Java knows which constructor to use
cg = CellGraph(regions, centroid_list, vertex_list, int32(0), int32(0), int32(6));  % 6 forces cells to be hexagons
verts = cg.vertexCoords;

%% build the parameter cell array
% fixed cells have at most 2 neighbors
fixed_cells_logical = cg.numOfCellsNeighboringVertex <= 2;
p.fixed_cells = find(fixed_cells_logical);
p.not_fixed_cells = find(~fixed_cells_logical);
p.numV = size(verts, 1);
p.preferred_distances = INITIAL_EQM_LENGTH_FACTOR * Misc.distMatrix(verts);
%This is where you can change the spring constant
p.spring_constants = SPRING_CONSTANT_INITIAL*ones(size(p.preferred_distances));

cgArray = javaArray('CellGraph', STEPS);
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
cg.setActiveCellsAuto(PASSIVE_LAYER_THICKNESS);

initial_preferred_distances = p.preferred_distances;
initial_cg = CellGraph(cg);

% so far connectivity doesn't change each round
% but with the stress reponse fibers it will
switch CONNECTIVITY
    case 'purse string'
        p.connectivity = cg.connectivityMatrixVertex;
    case 'network'
        p.connectivity = cg.connectivityMatrixVertexCrossConnections;
    case 'purse string and network'
        p.connectivity = cg.connectivityMatrixVertex + ...
        cg.connectivityMatrixVertexCrossConnections;        
    otherwise
        disp('No valid CONNECTIVITY selected');
end

%% main loop
disp('Beginning constriction...');
for j = 1:STEPS  
    % give the vertices a little kick
    verts = verts + JITTERING_STD*randn(size(verts));
    
    %% set active cells
    switch STRATEGY
        case 'synchronous'
            % keep all active Cells   
        case 'random asynchronous'
            cg.setActiveCells(cg.randomActiveCells(FRACTION_OF_ACTIVE_CELLS_TO_CONSTRICT));
        case 'neighbor asynchronous'
            if j == 1
                cg.setActiveCells(cg.randomActiveCells(FRACTION_OF_ACTIVE_CELLS_TO_CONSTRICT));
            else
                cg.setActiveCells(cg.activeNeighborsCells(cg.activeCells, FRACTION_OF_ACTIVE_CELLS_TO_CONSTRICT));
            end
        otherwise
            disp('No valid STRATEGY selected!');
    end

    %% make actconn, the connectivity matrix of active fibers
    if isempty(cg.activeCells)
        actconn = zeros(size(p.connectivity));
    else
        switch CONNECTIVITY
            case 'purse string'
                actconn = cg.connectivityMatrixVertexActiveCells;
            case'network'
                actconn = cg.connectivityMatrixVertexCrossConnectionsActiveCells;
            case 'purse string and network'
                actconn = cg.connectivityMatrixVertexActiveCells + ...
                cg.connectivityMatrixVertexCrossConnectionsActiveCells;
            otherwise
                disp('No valid CONNECTIVITY selected');
        end
    end
    actconn = logical(actconn);
    
    %% make preferred_distances

    % choose cells to constrict out of the active cells
    p.preferred_distances(actconn) = ...
        p.preferred_distances(actconn) * FRACTIONAL_NEW_EQUILIBRIUM_LENGTH;
    
    p.spring_constants(actconn) = ...
        p.spring_constants(actconn) * FRACTIONAL_NEW_SPRING_CONSTANT;

    p.initial_verts = verts;
    
    %% optimization routine
    % use the anonymous function @(x) calc_energy(x,p) to pass parameters...
    options = optimset('GradObj', 'on', 'DerivativeCheck', 'off', ...
        'LargeScale', 'off', 'Hessian', 'off');

    % remove the fixed cells from "verts"
    verts(p.fixed_cells, :) = [];
    % change the shape of verts to be a column vector (note: y before x)
    verts = [verts(:,1); verts(:,2)];
    [output_verts, E, exitflag, output, grad_true] = fminunc(@(x) calculate_energy(x, p), verts, options);
    output_verts = reshape(output_verts, length(output_verts)/2, 2);  % change it back
    % put the original vertices back in
    verts = p.initial_verts;
    verts(p.not_fixed_cells, :) = output_verts;
    
    % update cg with new vertex coordinates
    cg.updateVertexCoords(verts);
        
    %% reorganization phase
    switch ELASTICITY
        case 'reorganization'
            distmatrix = Misc.distMatrix(verts);
            p.preferred_distances = distmatrix;

%             switch CONNECTIVITY
%                 case 'purse string'
%                     reorconn = cg.connectivityMatrixVertex(cg.randomCells(FRACTION_OF_ACTIVE_CELLS_TO_CONSTRICT));
%                 case'network'
%                     reorconn = cg.connectivityMatrixVertexCrossConnections(cg.randomCells(FRACTION_OF_ACTIVE_CELLS_TO_CONSTRICT));
%                 case 'purse string and network'
%                     reorconn = cg.connectivityMatrixVertex(cg.randomCells(FRACTION_OF_ACTIVE_CELLS_TO_CONSTRICT)) + ...
%                         cg.connectivityMatrixVertexCrossConnections(cg.randomCells(FRACTION_OF_ACTIVE_CELLS_TO_CONSTRICT));
%                 otherwise
%                     disp('No valid CONNECTIVITY selected');
%             end
%             reorconn = logical(reorconn);
%             p.preferred_distances(reorconn) = distmatrix(reorconn);

        case 'elastic'
            % keep preferred_distances
        otherwise
            disp('No valid ELASTICITY selected!');
    end
    
    
    %% store the new CellGraph
    cgArray(j) = CellGraph(cg);
    % reset the active cells (if desired???)
%     cg.setActiveCellsAuto(PASSIVE_LAYER_THICKNESS)
    
%     %% calculate forces, anisotropy, cell areas, etc.
%     [bigvert biggrad stress] = calculate_forces(verts, p);
%     % note: grad gives x and y components of force- good for quiver!
%     
%     ellipse_properties = get_ellipse_properties(cg.drawActive);
%     major = ellipse_properties(:, 1);
%     minor = ellipse_properties(:, 2);
%     orientation = degtorad(ellipse_properties(:, 3));
%     
%     % anisotropy
%     r.average_anisotropy(j) = mean(major./minor);
%     
%     % x/y anisotropy
%     x0 = major .* minor ./ ... 
%         sqrt(major.^2.*sin(orientation).^2 + ...
%              minor.^2.*cos(orientation).^2);
%     y0 = major .* minor ./ ... 
%         sqrt(major.^2.*cos(orientation).^2 + ...
%              minor.^2.*sin(orientation).^2);
%     r.average_xy_anisotropy(j) = mean(x0 ./ y0);
%     
%     % constriction
%     areas = cg.activeCellAreas;
%     initial_areas = initial_cg.activeCellAreas;
%     r.average_constriction(j) = mean(areas ./ initial_areas);
%     
%     r.average_constriction_std(j) = ... 
%         std(areas ./ initial_areas)/sqrt(length(initial_areas));
%     
%     % rms stress
%     r.rms_stress(j) = sqrt(mean(biggrad(:).^2));
%     
%     % maximum stress
%     r.max_stress(j) = max(stress);
%     
%     % stress anisotropy
%     rms_AP_stress = sqrt(mean(biggrad(:, 2).^2));
%     rms_LR_stress = sqrt(mean(biggrad(:, 1).^2));
%     r.stress_anisotropy(j) = rms_AP_stress / rms_LR_stress;
%     r.stress_anis(j) = mean(biggrad(:,2)./biggrad(:,1));
%     
%     %% store values
%     P{j} = p;
%     Force{j, 1} = bigvert;
%     Force{j, 2} = biggrad;
%     Force{j, 3} = stress;
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

%% plotting
drawing(:,:,3) = initial_cg.draw;
drawing(:,:,1) = initial_cg.drawVertices;
figure('Position', [100 460 600 300]); 
imshow(drawing);

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

    %% plot active cells in yellow
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
