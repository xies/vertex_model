classdef Tissue < handle
    % TISSUE
    % Subclassed to handle.
    % Copy using new_tissue = Tissue(old_tissue) to make shallow copies
    %
    % ---- Properties ----
    %   vert_coords - Nv x 2 array of vertex locations
    %   vertices % hashmap of Vertex objects keyed by ID
    %   interfaces % hashmap of Interface.m objects
    %   cells - container.Map hashmap of CellModels keyed by cellID
    %   connectivity - Nv x Nv adjacency matrix of how vertices are
    %          connected
    %   interVertDist - Nv x Nv distance matrix of connected vertices
    %
    %   parameters - see setParameters for details
    %   contraction_param - contractile parameters
    %
    %   Ys - tissue image size (for display and initiation only)
    %   Xs - tissue image size
    %   merge_threshold_in_px - merge all vertices closer than this threshold
    %   t - time stamp
    %   dir - output directory
    %
    % ---- Methods ----
    % 
    %   Tissue - constructor
    %   isValid - consistency checker
    % 
    %   --- Energy methods ---
    %       get_energy - calculate the pot. energy of current config
    %       get_force - calculates the velocity on each vertex
    %
    %   --- Measurements ---
    %       getAreas - return area as function of sim time
    %       getContractility - contractility
    %       getCentroidX - x-coordinates
    %       getCentroidY - y-coordinates
    %       getTime - simulation time
    %   
    %   --- Simulation housekeeping methods ---
    %       move_vts - make new configuration via updating vert_coords
    %       moveVertex - move single vertex and update tissue
    %       t1Transition - T1 transition
    %       setParameters - sets the parameters and connects vertices
    %       updateVertCoords - update vert_coords explicitly
    %       activateCell - set certain cells to be active
    %       deactivateCell  - set certain cells to be not active
    %       deactivateBorder - set the outermost cells to be not active
    %       setContractility - directly set contractility by cellID
    %       setContractilityModel - set contractility based on model
    %          function of where cells are in the tissue
    %       sort - sort a tissue array by time
    %
    %   --- ODE solver methods ---
    %       step - one step of solver
    %       solve - calls ODE solver
    %
    %   --- Vertex-vertex connectivity ---
    %       connect_interfaces - sets or updates the 'edges'/'interfaces'
    %           of the model upon which energy/force is calculated
    %       edgeConnected - if two vertices are connected by edge
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
    %   --- Cell handling ---
    %       getCells - returns CellModels by ID or return all
    %       getActiveCells - return CellModels that are active
    %       getInactiveCells - return CellModels that are active
    %       setPerimElasticity - set perimeter elasticity to specified
    %          cells
    %       setAreaElasticity - set area elasticity to specified cells
    %       setTargetArea - set target area to specified cells
    %
    %   --- Interface handling ---
    %       setLineTension - set line tension
    %       setLineAnisotropy - set anisotropic line tension
    %       getInterfaces - get Interfaces by ID or return all
    %
    %   --- Vertex handling ---
    %       merge_vertices_initializing - merges vertices closer than a
    %          certain threshold (use during constructing only)
    %       make_vertices - gets rid of vertices beyond image, as well
    %           as return which cell touches a Vertex (only when
    %           constructing)
    %       getVertices - get Vertex by ID or return all
    %   
    %   --- Visualize ---
    %       draw - draw single tissue configuration in image
    %       movie - returns a (Xs x Ys x T) image stack of tissue
    %              configuration as a function of time
    %
    % xies@mit.edu March 2015
    
    properties
        
        % -- current config --
        vert_coords % Nvx2 array of vertex locations
        vertices % hashmap of Vertex.m objects
        cells % hashmap of the CellModels contained by the tissue
        interfaces % hashmap of Interface
        connectivity % adjacency matrix
        interVertDist % dist matrix
        parameters % simulation parameters
        contractile_params % contraction model + params and temporal model + params
        energy % energy
        t % timestamp
		t1Time % List of T1 transition times
		t1List % list of T1 transition vertices
        
        % -- Global settings --
        Xs % tissue pixel size
        Ys
        merge_threshold_in_px
        dir % output directory
        
    end
    methods
        
        function tis = Tissue(regions,vert_coords,centroids,conn_opt,t)
            % Constructor for Tissue object
            %
            % USAGE
            %  (to create from scratch):
            %       tis = Tissue(regions, vert_coords, centroids,'purse string',t)
            %       tis = Tissue(regions, vert_coords, centroids,'purse string')
            %                       assumes t = 0
            %  (to copy object)
            %        Tissue(old_tissue);
            %
            % March 2015, xies@mit.edu
            if nargin > 0
                if nargin < 5
                    t = 0;
                end
                
                if nargin == 1
                    % If there is a single input value, we must only want
                    % to copy an old tissue object (useful because
                    % container.Map is a reference object, and copies will
                    % modify the original)
                    
                    % Copy the value objects
                    tis_old = regions;
                    tis.Xs = tis_old.Xs; tis.Ys = tis_old.Ys;
                    tis.merge_threshold_in_px = tis_old.merge_threshold_in_px;
                    tis.t = tis_old.t;
                    tis.t1List = tis_old.t1List;
                    tis.t1Time = tis_old.t1Time;
                    tis.vert_coords = tis_old.vert_coords;
                    tis.connectivity = tis_old.connectivity;
                    tis.interVertDist = tis_old.interVertDist;
                    tis.parameters = tis_old.parameters;
                    tis.contractile_params = tis_old.contractile_params;
                    
                    % Copy the reference object via concatenation w/ empty
                    % Map, creating a brand new Map!
                    % DO NOT REMOVE
                    tis.cells = [tis_old.cells; containers.Map()];
                    tis.vertices = [tis_old.vertices; containers.Map()];
                    tis.interfaces = [tis_old.interfaces; containers.Map()];
                    
                else % Else construct from scratch
                    
                    tis.t = t;
                    tis.Xs = size(regions,1); tis.Ys = size(regions,2);
                    
                    tis.vertices = containers.Map('KeyType','int32','ValueType','any');

                    % Merge vertices that are too close to each other
                    % @todo: This requires speedup
                    tis.merge_threshold_in_px = 3;
                    vert_coords = tis.merge_vertices_initializing( ...
                        vert_coords, tis.merge_threshold_in_px);
                    tis.vert_coords = vert_coords;
                    
                    % Instantiate valid vertices
                    tis.make_vertices(regions);
                    
                    % Initialize a Map object to hold CellModel
                    num_cells = max(unique(regions));
                    tis.cells = containers.Map('KeyType','int32','ValueType','any');
                    
                    % Instatiate valid cells
                    vIDlist = tis.vertices.keys; vIDlist = [vIDlist{:}];
                    for i = 1:num_cells
                        tis.cells(i) = ...
                            CellModel(int32(i), tis,...
                            vIDlist( ...
                            cellfun(@(x) any(x == i), {tis.getVertices.cellIDs})), ... % Vertex keys
                            centroids(i,:) );
                    end
                    
                    % Make interfaces
                    tis.connect_interfaces(conn_opt);
                    tis.interVertDist = squareform( pdist(tis.vert_coords) );
                    tis.t1List = [NaN NaN]';
                    tis.t1Time = NaN;
                    
                end
                
                if ~isValid(tis)
                    error('Tissue not self-consistent, exiting');
                end
                
            end
        end % Constructor
        
        % Consistency checker
        flag = isValid(tis)
        
        function flags = isempty(tisArr)
            flags = false(1,numel(tisArr));
            flags = flags | cellfun(@isempty,{tisArr.cells});
%                 flag = 1; return;
%             end
%             if tis.vertices.length == 0,
%                 flag = 1; return;
%             end
%             if tis.interfaces.length == 0,
%                 flag = 1; return;
%             end
        end
        
        % ------  Calculate energy, force ------
        F = get_force(tis)
        E = get_energy(tis)
        
        % ------ Measurements (area, contractility, etc) ---------
        function A = getArea(tisArray)
            % Returns an array (or vector) of the cell areas contained in
            % input tissue array.
            T = numel(tisArray);
            % @todo: assumes tisArray does not change number of cells over
            % time!!
            % Idea for imporving: use cellID to make a NaN-padded array?
            num_cells = tisArray(1).cells.length;
            A = nan( T,num_cells );
            for i = 1:T
                c = tisArray(i).getCells;
                c = c.sortByID;
                A(i,:) = [c.area];
            end
        end % getArea
        
        function C = getContractility(tisArray)
            % Returns an array (or vector) of the contractility contained in
            % input tissue array.
            T = numel(tisArray);
            % @todo: assumes tisArray does not change number of cells over
            % time!!
            % Idea for imporving: use cellID to make a NaN-padded array?
            num_cells = tisArray(1).cells.length;
            C = nan( T,num_cells );
            for i = 1:T
                c = tisArray(i).getCells;
                c = c.sortByID;
                C(i,:) = [c.contractility];
            end
        end % getContractility
        
        function CX = getCentroidX(tisArray)
            % Returns an array (or vector) of the centroid-x contained in
            % input tissue array.
            T = numel(tisArray);
            % @todo: assumes tisArray does not change number of cells over
            % time!!
            % Idea for imporving: use cellID to make a NaN-padded array?
            num_cells = tisArray(1).cells.length;
            CX = nan( T,num_cells );
            for i = 1:T
                c = tisArray(i).getCells;
                c = c.sortByID;
                cx = cat(1,c.centroid);
                CX(i,:) = cx(:,2);
            end
        end % getCentroidX
        
        function CY = getCentroidY(tisArray)
            % Returns an array (or vector) of the centroid-y contained in
            % input tissue array.
            T = numel(tisArray);
            % @todo: assumes tisArray does not change number of cells over
            % time!!
            % Idea for imporving: use cellID to make a NaN-padded array?
            num_cells = tisArray(1).cells.length;
            CY = nan( T,num_cells );
            for i = 1:T
                c = tisArray(i).getCells;
                c = c.sortByID;
                cy = cat(1,c.centroid);
                CY(i,:) = cy(:,1);
            end
        end % getCentroidY
        
        function VX = getVertexX(tisArray)
            T = numel(tisArray);
            % @todo: assumes tisArray does not change number of vts over
            % time!!
            % Idea for imporving: use cellID to make a NaN-padded array?
            num_vertices = tisArray(1).vertices.length;
            VX = nan( T,num_vertices );
            for i = 1:T
                v = tisArray(i).getVertices;
                v = v.sortByID;
                VX(i,:) = cat(1,v.x);
            end
        end % getVertexX
        
        function VY = getVertexY(tisArray)
            T = numel(tisArray);
            % @todo: assumes tisArray does not change number of vts over
            % time!!
            % Idea for imporving: use cellID to make a NaN-padded array?
            num_vertices = tisArray(1).vertices.length;
            VY = nan( T,num_vertices );
            for i = 1:T
                v = tisArray(i).getVertices;
                v = v.sortByID;
                VY(i,:) = cat(1,v.y);
            end
        end % getVertexY

        function A = getAnisotropy(tisArray)
            % Returns an array (or vector) of the centroid-y contained in
            % input tissue array.
            T = numel(tisArray);
            % @todo: assumes tisArray does not change number of cells over
            % time!!
            % Idea for imporving: use cellID to make a NaN-padded array?
            num_cells = tisArray(1).cells.length;
            A = nan( T,num_cells );
            for i = 1:T
                c = tisArray(i).getCells;
                c = c.sortByID;
                A(i,:) = cat(1,c.anisotropy);
            end
        end % getAnisotropy
        
        function T = getTime(tisArray)
            % Returns a vector of the simulation time contained in
            % input tissue array.
            T = [tisArray.t];
        end % getTime
        
        % ------ Simulation bookkeeping methods ---------
        
        move_vts(tis, new_vcoords, new_time, varargin)
        setParameters(tis,varargin)
        t1Transition(tis,vt);
        
        function updateVertCoords(tis)
            % If you move a vertex without moving vert_coords explicitly,
            % use this to make them consistent.
            v = tis.getVertices;
            vx = [v.x]; vy = [v.y];
            tis.vert_coords = cat(2,vx',vy');
        end
        
        function moveVertex(tis,vID,new_coord)
            % Move a single vertex
            tis.vertices(vID) = tis.vertices(vID).move(new_coord);
            tis.updateVertCoords;
        end
        
        % -- Activate cell (mark as "ventral") --
        activateCell( tis, cellIDs, varargin)
        deactivateCell(tis, cellIDs)
        deactivateBorder(tis,n)
        
        % -- set contractility value to active cells only --
        setContractility(tis,C,cIDs)
        setContractilityModel(tis, contractions, cIDs)
        
        function jitterVertices(tis, STD)
            %JITTER_VERTICES
            % Add a set amount of Gaussian jitter to vertex position.
            %
            % USAGE: tis.jitterVertices( STD )
            
            verts = tis.vert_coords;
            I = ~tis.parameters.fixed_verts;
            jitter = STD*randn([ numel(I(I)), 2]);
            verts( I,: ) = verts( I,: ) + jitter;
            
            tis.move_vts( verts , 'no_update' );
            
        end % jitterVertices
        
        function tisArr = sort(tisArr)
            %SORT
            % Sort an array of Tissue objects by time stamp. Non-unique
            % time stamps will throw an error.
            % 
            % USAGE: tisArr = tisArr.sort;
            if numel(tisArr) == 1, return; end
            T = [tisArr.t];
            if numel(unique(T)) ~= numel(T)
                error('Non unique time stamp');
            end
            [~,I] = sort(T,'ascend');
            tisArr = tisArr(I);
        end % sort
        
        % --- ODE solver functions ---
        dy = step(tis,t,y,OUT_DIR)
        tisArr = solve(tis,tspan,OUT_DIR)
        
        % ------ Vertex-vertex connectivity ------
        
        connect_interfaces(tis, opt)
        function flag = edgeConnected(~,va,vb)
            % If two vertices are connected by an Interface
            flag = numel( intersect(va.bIDs,vb.bIDs) ) > 0;
        end % edgeConnected
        
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
        
        neighbors = neighborsOfCell( tis, cellOI, order_n )
        neighbors = allNthOrderNeighbors( tis, cellOI, order_n )
        
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
            va = cell_a.vIDs;
            for i = 1:numel(va)
                % compare vID / vertices key values
                flag = any(va(i) == cell_b.vIDs);
                if flag, return; end
            end
        end % connected
        
        function M = corona_measurement(tis,cellIDs,measurement)
            if nargin < 3,
                measurement = cellIDs;
                cellIDs = tis.cells.keys; cellIDs = [cellIDs{:}];
            end
            M = cell(1,numel(cellIDs));
            for i = 1:numel(cellIDs)
                neighbors = tis.neighborsOfCell(cellIDs(i),1);
                M{i} = [neighbors.(measurement)];
            end
        end % corona_measurement
        
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
                cellID = varargin{1};
                num_cells = numel(cellID);
                cells(1:num_cells) = CellModel();
                for i = 1:num_cells
                    cells(i) = tis.cells( cellID(i) );
                end
            end
        end % getCells
        
        function cells = getActiveCells(tis,returnID)
            % Return all active cells in the tissue (or just the IDs)
            %
            % USAGE: actives = tissue.getActiveCells;
            %        aIDs = tissue.getActiveCells(1);
            cells = tis.getCells;
            cells = cells( [cells.isActive] > 0 );
            if nargin < 2, returnID = 0; end
            if returnID, cells = [cells.cellID]; end
        end % getActiveCells
        
        function cells = getInactiveCells(tis)
            % Return all non-active cells in the tissue
            %
            % USAGE: actives = tissue.getInactiveCells;
            cells = tis.getCells;
            cells = cells( [cells.isActive] == 0 );
        end % getActiveCells
        
        function cIDs = getCellsWithinRegion(tis, box)
            % Returns all cells whose centroid lies within a region
            % specified by BOX = [X0 Y0 Xf Yf] (inclusive)
            %
            % USAGE: cellIDs = tis.getCellsWithinRegion( box )
            
            cells = tis.getCells;
            ct = cat(1,cells.centroid);
            cellIDList = tis.cells.keys();
            box([1,3]) = box([1,3])*tis.Xs;
            box([2,4]) = box([2,4])*tis.Ys;
            
            I = ct(:,1) >= box(1) & ct(:,1) <= box(3);
            I = I & ct(:,2) >= box(2) & ct(:,2) <= box(4);
            
            cIDs = [cellIDList{I}];
            
        end % getCellsWithinRegion
        
        % -- set parameter values to specific CellModel objects --
        setTargetArea(tis,cIDs,alt_area)
        setAreaElasticity(tis,cIDs,alt_ka)
        setPerimElasticity(tis,cIDs,alt_kp)
        
        % ----- Interface handling ----
        
        % -- set parameter values to specific Interface objects --
        setLineTension(tis,bIDList,L)
        setLineAnisotropy(tis,bIDList,a)
        
        function e = getInterfaces(tis, bIDs)
            % Returns interfaces from the tissue with specified IDs
            % Else returns all.
            % 
            % USAGE: e = tis.getInterfaces();
            % USAGE: e = tis.getInterfaces(bIDs);
            
            if nargin < 2,
                bIDs = tis.interfaces.keys;
                bIDs = [bIDs{:}];
            end
            
            num_edges = numel(bIDs);
            e(1:num_edges) = Interface;
            for i = 1:num_edges
                e(i) = tis.interfaces( bIDs(i) );
            end
        end % getInterfaces
        
        function e = getInterfaceByvIDs(tis,vIDs)
            % Returns the edge that has the given vIDs
            e = tis.getInterfaces;
            allVIDs = cellfun(@sort,{e.vIDs},'UniformOutput',0);
            e = e(cellfun(@(x) all( x==sort(vIDs) ),allVIDs));
            if numel(e) > 1, keyboard; end
        end % getInterfaceByvIDs
        
        % ----- Vertex handling ------
        
        make_vertices(tis, regions);
        vcoords = merge_vertices_initializing(~,vcoords,merge_threshold_in_px);
        
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
        
        draw(tis,varargin)
        vid = movie(tissues,varargin)
        cellID = select_cell(tis,n)
        vID = select_vertices(tis,n)
        
    end % Public methods
    
    
end
