classdef Tissue
    % TISSUE
    %
    % ---- Properties ----
    %   centroids - Nc x 2 array of centroid locations
    %   vert_coords - Nv x 2 array of vertex locations
    %   vertices % Nv x1 array of Vertex.m objects
    %   cells - container.Map hashmap of CellModels keyed by cellID
    %   connectivity - Nv x Nv adjacency matrix of how vertices are
    %          connected
    %   interVertDist - Nv x Nv distance matrix of connected vertices
    %
    %   --- Parameters - simulation parameters (see setParameters)
    %       p.targetAreas - target cell area
    %       p.targetPerimeter - target edge lengths
    %       p.areaElasticity - bulk elastic constant for cell area
    %       p.fixed_verts - list of fixed vertices
    %
    %   Ys - tissue image size (for display and initiation only)
    %   Xs - tissue image size
    %   merge_threshold_in_px - merge all vertices closer than this threshold
    %   t - time stamp
    %
    % ---- Methods ----
    % 
    %   Tissue - constructor
    %   isValid - consistency checker
    % 
    %   --- Energy methods ---
    %       get_energy - calculate the pot. energy of current config
    %       get_velocity - calculates the velocity on each vertex
    %   
    %   --- Simulation methods ---
    %       evolve - make new configuration via updating vert_coords
    %       serParameters - sets the parameters and connects vertices
    %       activateCell - set certain cells to be active
    %       deactivateCell  - set certain cells to be not active
    %       deactivateBorder - set the outermost cells to be not active
    %
    %   --- Vertex-vertex connectivity ---
    %       connectVertices - sets the 'edges'/'interfaces' of the model
    %             upon which energy/force is calculated
    %
    %   --- Vertex-cell connectivity ---
    %       cellsContainingVertex - returns all cells that contain input
    %              vertex
    %
    %   --- Cell-cell connectivity ---
    %       connected - whether two cells share one vertex
    %       allNthOrderNeighbors - returns the filled-in "halo" of
    %              n-th order neighbors of input cell
    %       neighborsOfCell - returns just the ring of n-th order neighbors
    %              of input cell
    %       numCellNeighbors - returns # of immediate neighbors of input
    %              cell
    % 
    %   --- Cell/vertex handling ---
    %       getActivateCells - return CellModels that are active
    %       mergeVertices - merges vertices closer than a certain threshold
    %          (use during constructing only)
    %       validateVertices - gets rid of vertices beyond image, as well
    %           as return which cell touches a vertex (only when
    %           constructing)
    %   
    %   --- Visualize ---
    %       draw - draw single tissue configuration in image
    %       movie - returns a (Xs x Ys x T) image stack of tissue
    %              configuration as a function of time
    %
    % xies@mit.edu March 2015
    properties
        
        % -- current config --
        centroids % Ncx2 array of centroid locations
        vert_coords % Nvx2 array of vertex locations
        vertices % array of Vertex.m objects
        cells % hashmap of the CellModels contained by the tissue
        connectivity % adjacency matrix
        interVertDist % dist matrix
        parameters % simulation parameters
        
        Xs % tissue pixel size
        Ys
        merge_threshold_in_px
        t % timestamp
        
    end
    methods
        
        function tis = Tissue(regions,vert_coords,centroids,t)
            % Constructor for Tissue object
            %
            % USAGE
            %  (to create from scratch):
            %       tis = Tissue(regions, vert_coords, centroids ,t)
            %       tis = Tissue(regions, vert_coords, centroids )
            %                       assumes t = 0
            %  (to copy object)
            %       tis = Tissue(old_tissue);
            %
            % March 2015, xies@mit.edu
            if nargin > 0
                if nargin < 4
                    t = 0;
                end
                
                if nargin == 1
                    % If there is a single input value, we must only want
                    % to copy an old tissue object (useful because
                    % container.Map is a reference object, and copies will
                    % modify the original)
                    
                    % Copy the value objects
                    tis_old = regions;
                    tis.centroids = tis_old.centroids;
                    tis.Xs = tis_old.Xs; tis.Ys = tis_old.Ys;
                    tis.merge_threshold_in_px = tis_old.merge_threshold_in_px;
                    tis.t = tis_old.t;
                    tis.vert_coords = tis_old.vert_coords;
                    tis.connectivity = tis_old.connectivity;
                    tis.interVertDist = tis_old.interVertDist;
                    tis.parameters = tis_old.parameters;
                    
                    % Copy the reference object via concatenation w/ empty
                    % Map, creating a brand new Map!
                    % DO NOT REMOVE
                    tis.cells = [tis_old.cells; containers.Map()];
                    tis.vertices = [tis_old.vertices; containers.Map()];
                    
                else % Else construct from scratch
                    
                    tis.t = t;
                    tis.centroids = centroids;
                    tis.Xs = size(regions,1); tis.Ys = size(regions,2);
                    
                    tis.vertices = containers.Map('KeyType','int32','ValueType','any');

                    % Merge vertices that are too close to each other
                    % @todo: This requires speedup
                    tis.merge_threshold_in_px = 6;
                    vert_coords = tis.merge_vertices( ...
                        vert_coords, tis.merge_threshold_in_px);
                    tis.vert_coords = vert_coords;
                    
                    % Instantiate valid vertices
                    tis = tis.make_vertices(regions);
                    
                    % Initialize a Map object to hold CellModel
                    num_cells = max(unique(regions));
                    tis.cells = containers.Map('KeyType','int32','ValueType','any');
                    % Instatiate valid cells
                    vIDlist = tis.vertices.keys;
                    vIDlist = [vIDlist{:}];
                    for i = 1:num_cells
                        tis.cells(i) = ...
                            CellModel(int32(i), tis,...
                            vIDlist( ...
                            cellfun(@(x) any(x == i), {tis.getVertices.cellIDs})), ... % Vertex keys
                            centroids(i,:) );
                    end
                end
                
                if ~isValid(tis)
                    error('Tissue not self-consistent, exiting');
                end
                
            end
        end % Constructor
        
        function flag = isValid(tis)
            % Run self-consistency tests on current tissue state
            % Checks the following:
            %   1) vertices and vert_coords match and are ordered correctly
            %   2) vertices and the unique set of cell-owned vertices also
            %      match (not ordered)
            %   3) dimensions of connectivity and interVertDist are
            %      consistent with vert_coords
            
            flag = 1;
            % vertices and vert_coords match
            vx = [tis.getVertices.x]; vy = [tis.getVertices.y];
            flag = flag & all(all(cat(2,vx',vy') == tis.vert_coords));
            
            % vertices and the set of cell vertices match
            cells = tis.getCells;
            vt = tis.getVertices( unique([cells.vIDs]) );
            vx = unique([vt.x]); vy = unique([vt.y]);
            flag = flag && all(vx == unique( [tis.getVertices.x] ) );
            flag = flag && all(vy == unique( [tis.getVertices.y] ) );
            
            if ~isempty( tis.connectivity)
                flag = flag && tis.vertices.length == mean(size(tis.connectivity));
                flag = flag && tis.vertices.length == mean(size(tis.interVertDist));
            end
            
        end
        
        % ------  Calculate energy, force, velocity ------
        
        function E = get_energy(tis)
            % GET_ENERGY Returns the current energy of the system
            % 
            % USAGE: E = get_energy(tis)
            %
            % Right now implements area elasticity, parameter elasticity,
            % and active contractility.
            
            p = tis.parameters; % get parameters
            
            dist = tis.interVertDist .* tis.connectivity;
            lineTensionTerm = nansum(nansum( triu(dist) ) * p.lineTension);
            
            current_areas = [tis.getCells.area];
            areaElasticTerm = nansum( p.areaElasticity.* ...
                ([tis.getCells.area] - p.targetAreas).^2 );
            perimElasticTerm = nansum( p.perimElasticity.* ...
                ([tis.getCells.perimeter] - p.targetPerimeters).*2 );
            
            C = [tis.getCells.contractility];
            activeContractionTerm = nansum( C .* current_areas );
            
            E = areaElasticTerm + perimElasticTerm + ...
                activeContractionTerm + lineTensionTerm;

        end % get_energy
        
        function V = get_force(tis)
            % GET_FORCE
            % Returns the forces as given by -grad(E) of current tissue
            % configuration.
            %
            % USAGE: V = tis.get_velocities;
            %
            % Currently implements, for all cells c, vertices i (and its neighbor j):
            %
            %   - area elasticity:
            %       Fa = -sum_c( k*(Ac-A0)*grad_i(Ac) )
            %   - perimeter elasticity
            %       Fp = -sum_c( -k*(Pc-P0)*grad_i(Pc)
            %   - line tension:
            %       Fl = -sum_i( sum_j(sij*grad_i(Dij)) )
            %
            % The graduent terms are implemented as Nagai (2000).
            
            % Grab data
            vcoords = tis.vert_coords;
            D = tis.interVertDist;
            conn = tis.connectivity;
            
            % Grap parameters
            gamma = tis.parameters.lineTension;
            kappa_a = tis.parameters.areaElasticity;
            kappa_p = tis.parameters.perimElasticity;
            
            % Initialize
            num_verts = size(vcoords,1);
            V = zeros( size(vcoords) );
            
            vIDList = tis.vertices.keys;
            vIDList = [vIDList{:}];
            % For now, loop through; vectorize later?
            for i = 1:num_verts
                
                % Find connected vertices
                vi = tis.vertices(vIDList(i));
                J = find(conn(i,:) == 1); % Idx of connected vertices
                num_neighbors = numel( vi.cellIDs );
                
                % Only give nonzero velocities for vertices w/ more than 2
                % cell neighbors (non-fixed cells)
                if num_neighbors > 2
                    
                    % Find conn cells (this is much faster)
                    neighbors(1:num_neighbors) = CellModel(); % initialize
                    for j = 1:num_neighbors
                        neighbors(j) = tis.cells(vi.cellIDs(j));
                    end
                    
                    % Intialize
                    line_tension_term = [0 0];
                    area_elastic_term = [0 0];
                    perim_elastic_term = [0 0];
                    active_contraction_term = [0 0];
                    
                    % Generate line tension term in direction of the
                    % inter-vertex vector
                    % Since loop is only a few elements long, it's faster
                    % than vectorized version
                    for j = 1:numel(J)
                        line_tension_term = line_tension_term - ...
                           gamma * (vcoords(i,:) - vcoords(J(j),:)) / D(i,J(j));
                    end
                    
%                     tis.draw('showVectors',{-line_tension_term,i},'showActive');
%                     keyboard
                    
                    % Go through all CELLS associated with current vertex,
                    % and calculate cell elasticity.
                    % NOTE that we can't just go through "edges" themselves
                    for this_cell = neighbors
                        
                        % Need to sort vertices counter-clockwise
                        sortedVt = tis.getVertices( this_cell.vIDs );
                        sortedVt = sortedVt.sort( this_cell.centroid );
                        
                        % Find this current vertex and its immediate
                        % neighboring vertices
                        I = find( vi == sortedVt );
                        I = wrap( [I-1 I I + 1] , numel(sortedVt)); % circularly index
                        r = [sortedVt(I).x ; sortedVt(I).y];
                        
                        R = [0 -1; 1 0]; % pi/2 rotation matrix
                        v = ( R*(r(:,1) - r(:,3)) )'; % grad(A)
                        
                        ua = r(:,1) - r(:,2); ub = r(:,3) - r(:,2);
                        u = ua / norm(ub) + ub /norm(ub); % Direction of grad(P)
                        
                        % Area elasticity
                        area_elastic_term = area_elastic_term ...
                            - 2 * (this_cell.area - tis.parameters.targetAreas) ...
                            * v * kappa_a;
                        
                        % Perimeter elasticity
                        perim_elastic_term = perim_elastic_term ...
                            + 2 * (this_cell.perimeter) * u' * kappa_p;
                        
                        % Active contraction -- same form as area
                        % elasticity but sets A0 = 0 and uses different
                        % coefficient
                        active_contraction_term = active_contraction_term ...
                            - 2 * this_cell.area * this_cell.contractility ...
                            * v;
                        
                        % DEBUGGING
%                         tis = tis.activateCell( this_cell.cellID );
%                         tis.draw('showVectors',{v,i},'showActive');
%                         tis = tis.deactivateCell( this_cell.cellID );
%                         hold on, sortedVt(I).draw;
%                         keyboard
                        
                    end
                    
                    V(i,:) = line_tension_term + area_elastic_term + ...
                        perim_elastic_term + active_contraction_term;
                    
%                     tis.draw('showVectors',{V(i,:),i},'showActive');
%                     drawnow

                end
                
                if any(any(isnan( V ))), keyboard; end
                
            end
            
            % Check that fixed vertices have not moved
            if any( V(tis.parameters.fixed_verts,:) > 0 )
                keyboard;
            end
            
        end % get_velocities
        
        % ------ Simulation methods ---------
        
        function tis = evolve( tis_old, new_vcoords, varargin)
            % EVOLVE
            % Updates and returns a new copy of the old tissue
            % configuration by moving all the vertex positions
            %   NOTA BENE: the old .cells container is COPIED and not
            %              direclty modified, since it's a reference object
            %
            % USAGE: new_tissue = old_tissue( new_vcoords );
            %        new_tissue = old_tissue( new_vcoords ,'no_update');
            
            if size(new_vcoords,1) ~= tis_old.vertices.length
                error('Size of new vertex list must match old vertices')
            end
            tis = Tissue(tis_old);
            tis.vert_coords = new_vcoords;
            vIDList = tis.vertices.keys; vIDList = [vIDList{:}];
            
            for i = 1:numel(vIDList)
                % Move vertices in CellModels
                v = tis.vertices(vIDList(i));
                cContainV = v.cellIDs;
                for c = cContainV'
                    tis.cells( c ) = ...
                        tis.cells(c).updateCell( tis );
                end
                
                % Move Vertex
                tis.vertices(vIDList(i)) = ...
                    tis.vertices(vIDList(i)).move(new_vcoords(i,:));
            end
            
            % Update distance maps
            tis.interVertDist = squareform(pdist(tis.vert_coords));
            
            % Consistency check
            if ~tis.isValid, keyboard; end
            
            % Advance time stamp by one
            if nargin > 2,
                if ~strcmpi(varargin{1},'no_update'); tis.t = tis.t + 1; end
            end
            
        end % evolve
        
        function tis = setParameters(tis,parameters)
            % Sets the simluation/evolution parameters
            %
            % USAGE: tissue = 
            %           tis.setParameters(p);
            % 
            % INPUT: tis - tissue
            %        p.targetArea - target area
            %        p.lineTension - line tension
            %        p.areaElasticity - area elasticity
            %        p.connect_opt - connectivity option ('purse string')
            %
            % @todo: Figure out how to error-handle bad inputs
            
            tis.parameters = parameters;
            tis.parameters.fixed_verts = cellfun( ...
                @numel,{tis.getVertices.cellIDs}) < 3;
            conn = tis.adjMatrix(parameters.conn_opt);
            tis.connectivity = conn;
            tis.interVertDist = squareform( pdist(tis.vert_coords) );
            
        end % setParameters
        
        function tis = activateCell( tis, cellIDs, varargin)
            % activateCell
            % Set specified cells (IDs) to "active = 1"
            % Usage: tissue = activateCell(tissue); Activates all cells
            %        tissue = activateCell(tissue, [1 2 3]); Activates the
            %              specified subset of cells
            %        tissue = activateCell('random',fraction); Activates
            %              random subset of cells
            % AVAILABLE STRATEGIES:
            %        'random' - randomly activate up to the given fraction
            
            % If no arguments are given, then activate all cells
            if nargin < 2,
                cellIDs = tis.cells.keys;
                cellIDs = [cellIDs{:}];
            % If more than 2 arguments are given, then need to set the
            % cellIDs
            elseif nargin > 2
                strategy = cellIDs;
                switch strategy
                    case 'random'
                        % Randomly activate up to specified fraction
                        fraction = varargin{1};
                        if fraction > 1, fraction = 1; warning('fraction > 1'); end
                        num_cells = tis.cells.length;
                        ones2Activate = false(1,num_cells);
                        ones2Activate( 1:round(fraction*num_cells) ) = true;
                        ones2Activate = ones2Activate( randperm(num_cells) );
                        
                        cellIDs = tis.cells.keys();
                        cellIDs = [cellIDs{ones2Activate}];
                    otherwise
                end
            end
            % Activate specified cells
            for i = 1:numel(cellIDs)
                tis.cells( cellIDs(i) ) = ...
                    tis.cells( cellIDs(i) ).activateCell;
            end
        end % activateCell
        
        function tis = deactivateCell(tis, cellIDs)
            % Deactivates cell(s)
            % USAGE: tis = tis.deactivateCell; (default = all cells)
            %        tis = tis.deactivateCell(IDs);
            if nargin < 2 % If no cells specified, deactivate all cells
                cellIDs = tis.cells.keys;
                cellIDs = [cellIDs{:}];
            end
            for i = 1:numel(cellIDs)
                tis.cells( cellIDs(i) ) = ...
                    tis.cells( cellIDs(i) ).deactivateCell;
            end
        end % deactivateCell
        
        function tis = deactivateBorder(tis,n)
            % deactivateBorder:
            % Deactivates the border (up to n order) cells
            % USAGE: tis = tis.deactivateBorder(n)
            %
            % @todo: implement orders n > 1
            
            % By default only deactivate first layer
            if nargin < 2, n = 1; end
            
            if n>1, error('Haven''t implemented higher order :('); end
            
            cellIDList = tis.cells.keys; cellIDList = [cellIDList{:}];
            for ID = cellIDList
                if tis.numCellNeighbors( tis.cells(ID) ) < 6
                    tis = tis.deactivateCell( ID );
                end
            end
            
        end % deactivateBorder
        
        function tis = setContractility(tis,C)
            % Directly sets the active contractility coefficient
            % Requires all cells to have its own specified contractility
            % 
            % USAGE: tis = tis.setContractility( C );
            % INPUT: tis - tissue
            %        C - (Nc x 1) vector of contractility values
            
            if tis.cells.length ~= numel( C )
                error('Number of contractility coeff and number of cells don''t match');
            end
            
            cellIDList = tis.cells.keys;
            for i = 1:tis.cells.length
                tis.cells( cellIDList{i} ) = ...
                    tis.cells( cellIDList{i} ).setContractility(C(i));
                    % @todo: need to figure out edge-contractiltiy and how
                    % to inherit it from a cell
                
            end
            
        end
        
        function tis = jitterVertices(tis, STD)
            % Add a set amount of Gaussian jitter to vertex position.
            %
            % USAGE: tis = tis.jitterVertices( STD )
            
            verts = tis.vert_coords;
            I = ~tis.parameters.fixed_verts;
            jitter = STD*randn([ numel(I(I)), 2]);
            verts( I,: ) = verts( I,: ) + jitter;
            
            tis = tis.evolve( verts , 'no_update' );
            
        end
        
        % ------ Verted-vertex connectivity ------
        
        function conn = adjMatrix(tis, opt, cells)
            % Sets the vertex-vertex connectivity matrix according to the
            % specified model configurations
            % 
            % USAGE: 
            % tis = tis.connectVertices( opt )
            % tis = tis.connectVertices( opt , cells )
            %      (only connect a subset of cells)
            % 
            % INPUT: tis - the tissue to be connected
            %        opt - 'purse string' / 'apical' / 'both'
            %        cells - cellIDs (@todo: not implemented fully)
            % 
            % @todo: implement 'apical' and 'both'
            
            vt = tis.getVertices;
            num_vertices = numel(vt);
            if nargin < 3
                cells = tis.getCells;
            end
            
            switch opt
                case 'purse string'
                    % Connect the 'interfaces' of cells only
                    % @todo: SLOWWW!!!
                    conn = zeros(num_vertices);
                    
                    for i = 1:num_vertices
                        for this_cell = cells
                            connectedVertIDs = this_cell.getConnectedVertices( vt(i) );
                            I = ismember([vt.ID], connectedVertIDs);
                            conn(i,I) = 1;
                        end
                    end
                    
                otherwise
                    error('Unrecognized vertex connection option.')
            end
            conn( logical(eye(num_vertices)) ) = 0;
            
        end % connectVertices
        
        % ------ Cell-Vertex connectivity -----
        
        function cellsThatTouch = cellsContainingVertex(tis,vert)
            % Return a list of CellModel that contains the current Vertex
            %
            % Usage: touchingCells = tis.cellsContainingVertex( vert )
            
            if numel(vert) > 1, error('Can''t do more than 1 vertex.'); end
            N = numel(vert.cellIDs);
            cellsThatTouch(1:N) = CellModel();
            for i = 1:N
                cellsThatTouch(i) = tis.cells( vert.cellIDs(i) );
            end
            
            if isempty(cellsThatTouch), keyboard; end
        end % cellsContainingVertex
        
        % ----- Cell-cell connectivity ----
        
        function neighbors = neighborsOfCell( tis, cellOI, order_n )
            % Return a annulus of cells around cellOI, at order_n-th layer
            % USAGE: neighbors = neighborsOfCell( tis, cellOI, order_n )
            
            % If input is only a cellID instead of object itself, then fetch cell
            if isnumeric(cellOI), cellOI = tis.cells(cellOI); end
            
            tmpCellSet = java.util.HashSet;
            % Add the n-th neighbors
            nth_neighb = tis.allNthOrderNeighbors(cellOI, order_n);
            for c = nth_neighb
                tmpCellSet.add( c.cellID );
            end
            % Delete the n-1-th neighbors
            nminus1_neighbor = tis.allNthOrderNeighbors(cellOI, order_n-1);
            for c = nminus1_neighbor
                tmpCellSet.remove( c.cellID );
            end
            
            % Grab neigbor cells from the Map in Tissue
            tmpCellSet = tmpCellSet.toArray;
            neighbors(1:numel(tmpCellSet)) = CellModel;
            for i = 1:numel(tmpCellSet)
                neighbors(i) = tis.cells(tmpCellSet(i));
            end
            
        end %neighborsOfCell
        
        function neighbors = allNthOrderNeighbors( tis, cellOI, order_n )
            % allNthOrderNeighbors - Find the "halo" of cells around cellOI
            % until we've filled it to ORDER_N layers.
            %
            % USAGE: neighbors = allNthOrderNeighbors( tis, cellOI, orderN )
            %
            % OUTPUT: neighbors - cellModel array
            %
            % Not very fast -- @todo: need HashSet to work with matlab custom
            % classes, can't figure out yet.
            
            % If input is only a cellID instead of object itself, then fetch cell
            if isnumeric(cellOI), cellOI = tis.cells(cellOI); end
            
            tmpCellSet = java.util.HashSet; % HashSet constant add time
            % Add the central cell; can only use primitive data types with
            % imported HashSet :(
            cellIDList = tis.cells.keys();
            tmpCellSet.add(cellOI.cellID);
            for i = 1:order_n
                
                % Grab all the current i-th order filled-in corona
                tmpCells = tmpCellSet.toArray;
                for d = 1:numel(tmpCells)
                    % Iterate through all cells to find all connected cells
                    thisCell = tis.cells( tmpCells(d) );
                    %                     if i == 2, keyboard; end
                    for c = 1:tis.cells.length
                        if tis.connected( tis.cells(cellIDList{c}), thisCell )
                            tmpCellSet.add( tis.cells(cellIDList{c}).cellID );
                        end
                    end
                end
                
            end
            
            % Grab neigbor cells from the Map in Tissue
            tmpCellSet = tmpCellSet.toArray;
            neighbors(1:numel(tmpCellSet)) = CellModel;
            for i = 1:numel(tmpCellSet)
                neighbors(i) = tis.cells(tmpCellSet(i));
            end
            
        end % allNthOrderNeighbors
        
        function numNeigh = numCellNeighbors(tis, cellA)
            % Return the number of neighboring cells to input cell
            % USAGE: numNeighbor = tis.numCellNeighbors(input_cell)
            numNeigh = numel( tis.neighborsOfCell( cellA,1 ) );
        end % numCellNeighbors
        
        function flag = connected(~,cell_a, cell_b)
            % Basic function where if two cells share one vertex,
            % they're connected
            % 
            % Usage: flag = tis.connected( cells_a, cell_b)
            
            flag = 0;
            va = cell_a.vertices;
            for i = 1:numel(va)
                flag = any(va(i) == cell_b.vertices);
                if flag, return; end
            end
        end % connected
        
        % ----- Cell handling and measurements ------
        
        function cells = getCells( tis, varargin )
            % Returns cells from tissue as an array.
            % 
            % USAGE:
            %    cells = tis.getCells(); % returns all cells
            %    cells = tis.getCells(cellID); returns subset
            
            % Get all cells
            if nargin == 1
                cells = tis.cells.values;
                cells = [cells{:}];
            else
                cellID = varargin;
                num_cells = numel(cellID);
                cells(1:num_cells) = CellModel();
                for i = 1:num_cells
                    cells(i) = tis.cells( cellID(i) );
                end
            end
        end % getCells
        
        function cells = getActiveCells(tis)
            % Return all active cells in the tissue
            %
            % USAGE: actives = tissue.getActiveCells;
            cells = tis.getCells;
            cells = cells( [cells.isActive] > 0 );
        end % getActiveCells
        
        % ----- Vertex handling ------
        
        function tis = make_vertices(tis, regions)
            % Checks if vertices are not beyond image border, touching at 
            % least one cell, and return which cells a vertex is touching
            % Use only from Constructor!
            
			% Check that tis.vertices is empty -- or else this is not
			% called from the Constructor.
			if tis.vertices.length > 0
				error('Don''t call this function outside of construtor');
			end

			vert_coords = tis.vert_coords;
            num_vertices = size( vert_coords,1 );
            for i = num_vertices:-1:1
                
                x = vert_coords(i,1); y = vert_coords(i,2);
            	
				% Check if vx is beyond boundary, if so, chuck it and
            	% move on
                if x < 1 || x > tis.Xs || y < 1 || y > tis.Ys
                    vert_coords(i,:) = [];
                    continue;
                end
                
                % Extract 8 connected neighbors of this vertex
                conn_pixels = regions(round(x)-1:round(x)+1, ...
                    round(y)-1:round(y)+1);
                neighbor_cells = unique(conn_pixels(conn_pixels > 0));
                
                % Check that vertex has at least 1 connected cell, if no,
                % then chuck it
                if isempty(neighbor_cells)
                    vert_coords(i,:) = [];
                else
                    % If is not empty- then instantiate a Vertex
                    % and add it to the Map (handle object so should
                    % be modified in all spaces)
                    tis.vertices( int32(i) ) = ...
                        Vertex( int32(i), vert_coords(i,1), vert_coords(i,2), ...
                        neighbor_cells);
                end
            end
            
            tis.vert_coords = vert_coords;
            
        end % make_vertices
        
        function vcoords = merge_vertices(~,vcoords, ...
                merge_threshold_in_px)
            % Merge vertices that are closer than the specified threshold
            %
            % USAGE:
            % vcoords = tis.merge_vertices( vcoords,threshold)
            %
            % Use only from constructor
            
            num_vertices = size(vcoords);
            % Merge vertices which are super close to each other
            vertDist = squareform(pdist(vcoords));
            vertDist( logical(eye(num_vertices)) ) = NaN;
            
			% Loop merge function until no pairs of vertices are closer than
			% the threshold
            while any(any(vertDist <= merge_threshold_in_px))
                
                % Find a set of vertices to merge
                [I,J] = find(vertDist <= merge_threshold_in_px,1,'first');
                tobeMergedInd = [I,find(vertDist(I,:) <= merge_threshold_in_px)];
				mergedV = [mean( vcoords(tobeMergedInd,:) )];
                
                % Delete old ummerged vertices
                vcoords(tobeMergedInd,:) = [];
                % Add new merged vertex
                vcoords = cat(1,vcoords,mergedV );
                
                % update the distance maps
                num_vertices = size(vcoords,1);
                vertDist = squareform(pdist(vcoords));
                vertDist( logical(eye(num_vertices)) ) = NaN;
                
            end
        end % merge_vertices
        
        function verts = getVertices( tis , IDs )
            % Returns all or specified subset of vertices as array
            % USAGE: vts = tis.getVertices
            % 	     vts = tis.getVertices( IDs )
            
            % Return all Vertex objects as an array
            if nargin < 2
                verts = tis.vertices.values;
                verts = [verts{:}];
            else
                verts(1:numel(IDs)) = Vertex();
                for i = 1:numel(IDs)
                    verts(i) = tis.vertices( IDs(i) );
                end
            end
        end % getVertices
        
        %------ Visualization -----
        
        function I = draw(tis,varargin)
            % Draws a single tissue in binary image. Can return just the
            % outlines of cells, or also shade-in the active cells.
            %
            % USAGE: I = draw(tis);
            %        I = draw(tis,'showActive');
            %           Highlights currentlly "active" cells
            %        I = draw(tis,'showVectors',V);
            %           Shows vector field V over all vertices
            %        I = draw(tis,'showVectors',{v,ID});
            %           Shows only single vector v at vertex == ID
            %
            % Cannot handle more than one tissue.
            
            if numel(tis) > 1, error('Can only handle single tissue; use tis.movie() to show movie.'); end
            
            I = zeros(tis.Xs,tis.Ys);
            cellIDList = tis.cells.keys();
            for i = 1:numel(cellIDList)
                I = I + tis.cells(cellIDList{i}).draw(tis);
                I = logical(I);
            end
            
            I = double(I) * 255;
            % Highlight active cells
            if any(strcmpi(varargin, 'showActive'))
                % Show active cells as filled-ins
                M = zeros(tis.Xs,tis.Ys);
                Acells = tis.getActiveCells;
                for i = 1:numel(Acells)
                    M = M + Acells(i).drawMask(tis);
                end
                M = M * 50;
                I = I + M;
            end
            
            imagesc(I), axis equal;
            
            % Display vector field
            ind = find( strcmpi(varargin,'showVectors') );
            if ~isempty(ind)
                if numel(varargin) < ind + 1, error('Need vector to draw'); end
                V = varargin{ind+1};
                % If input is not a cell object, then display all vectors
                if ~iscell(V)
                    hold on;
                    quiver(tis.vert_coords(:,2),tis.vert_coords(:,1), ...
                        V(:,2),V(:,1),0,'w-');
                else
                    % If it's a cell obj, then only the given vertex will
                    % have a vector over it
                    v = V{1};
                    ID = V{2};
                    if numel(ID) ~= size(v,1);
                        error('# of vectors should equal # of origins')
                    end
                    hold on
                    quiver(tis.vert_coords(ID,2),tis.vert_coords(ID,1), ...
                        v(:,2),v(:,1),0,'w-');
                end
            end
            
        end
        
        function F = movie(tissues,varargin)
            % Make a movie of tissue evolving. Can return just the
            % outlines of cells, or also shade-in the active cells.
            %
            % USAGE: I = movie(tisues);
            %        I = movie(tisues,'showActive');
            
            num_frames = numel(tissues);
            if nargin == 1, opt = 'none';
            else
                opt = varargin{1};
            end
            
            F = zeros(tissues(1).Xs, tissues(1).Ys, num_frames);
            
            for f = 1:num_frames
                F(:,:,f) = tissues(f).draw(opt);
            end
            
        end
        
    end % Methods
    
    
end
