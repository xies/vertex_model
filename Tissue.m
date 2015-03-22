classdef Tissue
    
    properties
        
        centroids % Ncx2 array of centroid locations
        vert_coords % Nvx2 array of vertex locations
        vertices % array of Vertex.m objects
        connectivity % Nv x Nv matrix of how vertices are connected
        cells % hashmap of the CellModels contained by the tissue
        
        Xs % tissue pixel size
        Ys
        merge_threshold_in_px
        t
        
    end
    methods
        
        function tis = Tissue(regions,vert_coords,centroids,t)
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
                    tis.vertices = tis_old.vertices;
                    tis.vert_coords = tis_old.vert_coords;
                    
                    % Copy the reference object via concatenation w/ empty
                    % Map, creating a brand new Map!
                    % DO NOT REMOVE
                    tis.cells = [tis_old.cells; containers.Map()];
                    
                else % Else construct from scratch
                    
                    tis.t = t;
                    tis.centroids = centroids;
                    tis.Xs = size(regions,1); tis.Ys = size(regions,2);
                    num_cells = max(unique(regions));
                    tis.cells = containers.Map('KeyType','int32','ValueType','any');
                    
                    % Get vertices that are not on border
                    [vert_coords,~] = tis.validate_vertices(regions,vert_coords);
                    num_vertices = size(vert_coords,1);
                    vertices(1:num_vertices) = Vertex; % preallocate empties
                    for i = 1:num_vertices
                        vertices(i) = Vertex(vert_coords(i,1),vert_coords(i,2));
                    end
                    
                    
                    % Merge vertices that are too close to each other
                    % @todo: This requires speedup
                    tic
                    tis.merge_threshold_in_px = 3;
                    [vertices,vert_coords] = tis.merge_vertices(vert_coords,vertices,...
                        tis.merge_threshold_in_px);
                    tis.vert_coords = vert_coords;
                    tis.vertices = vertices;
                    toc
                    
                    % Get cell-ownership of vertices via 8-connected neighbors of
                    % vertices and REGIONS map
                    [tis.vert_coords,vx2Cell] = tis.validate_vertices(regions);
                    
                    % Instatiate valid cells
                    for i = 1:num_cells
                        tis.cells(int32(i)) = ...
                            CellModel(int32(i), tis,...
                            vertices( cellfun(@(x) any(x == i),vx2Cell) ), ...
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
            flag = 1;
            % vertices and vert_coords match
            vx = [tis.vertices.x]; vy = [tis.vertices.y];
            flag = flag & all(all(cat(2,vx',vy') == tis.vert_coords));
            
            % vertices and the set of cell vertices match
            cells = tis.cells.values;
            cells = [cells{:}];
            vt = [cells.vertices];
            vx = unique([vt.x]); vy = unique([vt.y]);
            flag = flag && all(vx == unique( [tis.vertices.x] ) );
            flag = flag && all(vy == unique( [tis.vertices.y] ) );
            
        end
        
        % ------ Simulation methods ---------
        
        function tis = evolve( tis_old, new_vcoords)
            % NOTA BENE: Since tis.cells is a reference object, it's
            % dynamically updated, so be careful modifying its values
            
            if size(new_vcoords,1) ~= numel(tis_old.vertices)
                error('Size of new vertex list must match old vertices')
            end
            tis = Tissue(tis_old);
            tis.vert_coords = new_vcoords;
            
            for i = 1:numel(tis.vertices)
                
                % Move vertices in CellModels
                v = tis.vertices(i);
                cContainV = [tis.cellsContainingVertex(v).cellID];
                for c = cContainV
                    tis.cells( c ) = ...
                        tis.cells( c ).moveVertex(v,new_vcoords(i,:));
                end
                
                % Move Vertex
                tis.vertices(i) = tis.vertices(i).move(new_vcoords(i,:));
            end
            % Advance time stamp by one
            tis.t = tis.t + 1;
        end
        
        function tis = activateCell( tis, cellIDs, varargin)
            % Set specified cells (IDs) to "active = 1"
            % Usage: tissue = activateCell(tissue, [1 2 3]);
            %        tissue = activateCell(tissue); Activates all cells
            %        tissue = activateCell('random',fraction); Activates
            %               random subset of cells
            % AVAILABLE STRATEGIES:
            %        'random' - randomly activate
            %
            
            if nargin < 2,
                cellIDs = tis.cells.keys;
                cellIDs = [cellIDs{:}];
            elseif nargin > 2
                strategy = cellIDs;
                switch strategy
                    case 'random'
                        % Randomly activate up to specified fraction
                        fraction = varargin{1};
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
            % Deactivates the border (up to n order) cells
            % USAGE: tis = tis.deactivateBorder(1)
            
            % By default only deactivate first layer
            if nargin < 2, n = 1; end
            
            if n>1, error('Haven''t implemented higher order :('); end
            
            cellIDList = tis.cells.keys; cellIDList = [cellIDList{:}];
            for ID = cellIDList
                if tis.numCellNeighbors( tis.cells(ID) ) < 6
                    tis = tis.deactivateCell( ID );
                end
            end
            
        end
            
        
        % ------ Verted-vertex connectivity ------
        
        function tis = connectVertices(tis,opt, cells)
            %
            vt = tis.vertices;
            num_vertices = numel(vt);
            if nargin < 3
                cells = tis.cells.values;
                cells = [cells{:}];
            end
            
            switch opt
                case 'purse string'
                    % Connect the 'interfaces' of cells only
                    conn = zeros(num_vertices);
                    
                    for i = 1:num_vertices
                        for this_cell = cells
                            neighbors = this_cell.getConnectedVertices( vt(i) );
                            I = vt.ismember( neighbors );
                            conn(i,I) = 1;
                        end
                    end
                    
                otherwise
                    error('Unrecognized vertex connection.')
            end
            tis.connectivity = conn;
        end
        
        % ------ Cell-Vertex connectivity -----
        
        function cellsThatTouch = cellsContainingVertex(tis,vert)
            % Return a list of CellModel that contains the current Vertex
            cellsThatTouch = [];
            keylist = tis.cells.keys();
            for key = keylist
                c = tis.cells(key{:});
                if any([c.vertices] == vert)
                    cellsThatTouch = [cellsThatTouch c];
                end
            end
            if isempty(cellsThatTouch), keyboard; end
        end % cellsContainingVertex
        
        function num_touch = numCellTouchingVertices(tis)
            % Return a list of CellModel that contains the current Vertex
            num_vertices = numel(tis.vertices);
            num_touch = zeros(1,num_vertices);
            % Go through all vertices
            for v = 1:num_vertices
                vert = tis.vertices(v);
                % Go through all cells and find ones containing vert
                c = tis.cellsContainingVertex( vert );
                num_touch(v) = numel(c);
            end
            
        end %numCellTouchingVertices
        
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
            
        end %allNthOrderNeighbors
        
        function numNeigh = numCellNeighbors(tis, cellA)
            % Return the number of neighboring cells to input cell
            % USAGE: numNeighbor = tis.numCellNeighbors(input_cell)
            numNeigh = numel( tis.neighborsOfCell( cellA,1 ) );
        end
        
        function flag = connected(~,cells_a, cells_b)
            % Basic function where if two cells share one vertex,
            % they're connected
            flag = 0;
            va = cells_a.vertices;
            for i = 1:numel(va)
                flag = any(va(i) == cells_b.vertices);
                if flag, return; end
            end
        end % connected
        
        %         function flag = connecte
        
        % ----- Cell handling ------
        
        function cells = getActiveCells(tis)
            % Return all active cells in the tissue
            %
            % USAGE: actives = tissue.getActiveCells;
            cells = tis.cells.values;
            cells = [cells{:}];
            cells = cells( [cells.isActive] > 0 );
        end
        
        % ----- Vertex handling ------
        
        function [vert_coords,vx2Cell] = validate_vertices(tis, regions, vert_coords)
            
            if nargin < 3, vert_coords = tis.vert_coords; end
            
            % Check if vx is beyond boundary, if so, chuck it and
            % move on
            num_vertices = size(vert_coords,1);
            for i = num_vertices:-1:1
                
                x = vert_coords(i,1); y = vert_coords(i,2);
                if x < 1 || x > tis.Xs || y < 1 || y > tis.Ys
                    vert_coords(i,:) = [];
                    continue;
                end
                
                % Extract 8 connected neighbors of this vertex
                conn_pixels = regions(x-1:x+1, y-1:y+1);
                neighbor_cells = unique(conn_pixels(conn_pixels > 0));
                vx2Cell{i} = neighbor_cells;
                
                % Check that vertex has at least 1 connected cell, if no,
                % then chuck it
                if isempty(vx2Cell{i})
                    vert_coords(i,:) = [];
                    vx2Cell(i) = [];
                end
                
            end
        end % validate_vertices
        
        function [vertices,vcoords] = merge_vertices(~,vcoords,vertices,...
                merge_threshold_in_px)
            
            num_vertices = size(vcoords);
            % Merge vertices which are super close to each other
            vertDist = squareform(pdist(vcoords));
            vertDist( logical(eye(num_vertices)) ) = NaN;
            while any(any(vertDist <= merge_threshold_in_px))
                
                % Find a set of vertices to merge
                [I,J] = find(vertDist <= merge_threshold_in_px,1,'first');
                tobeMergedInd = [I,find(vertDist(I,:) <= merge_threshold_in_px)];
                mergedV = vertices( tobeMergedInd ).merge;
                
                % Delete old ummerged vertices
                vertices(tobeMergedInd) = [];
                vcoords(tobeMergedInd,:) = [];
                % Add new merged vertex
                vertices = [vertices mergedV];
                vcoords = cat(1,vcoords, [mergedV.x mergedV.y] );
                
                % update the distance maps
                num_vertices = size(vcoords,1);
                vertDist = squareform(pdist(vcoords));
                vertDist( logical(eye(num_vertices)) ) = NaN;
                
            end
        end % merge_vertices
        
        %------ Visualization -----
        
        function I = draw(tis,opt)
            % Returns a binary image of the current tissue state
            if nargin < 2
                opt = 'none';
            end
            I = zeros(tis.Xs,tis.Ys);
            cellIDList = tis.cells.keys();
            for i = 1:numel(cellIDList)
                I = I + tis.cells(cellIDList{i}).draw;
                I = logical(I);
            end
                
            % Show active cells as filled-ins
            M = zeros(tis.Xs,tis.Ys);
            switch opt
                case 'none'
                    I = double(I) * 255;
                case 'showActive'
                    Acells = tis.getActiveCells;
                    for i = 1:numel(Acells)
                        M = M + Acells(i).drawMask;
                    end
                    M = M * 50;
                    I = double(I)*255 + M;
                otherwise
                    error('Unrecognized draw command/option')
            end
            
        end
        
        function F = movie(tissues,opt)
            % Make a movie of tissue evolving
            num_frames = numel(tissues);
            
            F = zeros(tissues(1).Xs, tissues(1).Ys, num_frames);
            
            for f = 1:num_frames
                F(:,:,f) = tissues(f).draw(opt);
            end
            
        end
        
    end % Methods
    
    
end
