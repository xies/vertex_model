classdef Tissue
    
    properties
        
        centroids % Nx2 array of centroid locations
        vert_coords % Nvx2 array of vertex locations
        vertices % array of Vertex.m objects
        regions % output of bwlabel, -1 = boundary, cells = 1,2,3...
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
                
                tis.t = t;
                tis.regions = regions;
                tis.centroids = centroids;
                tis.Xs = size(regions,1); tis.Ys = size(regions,2);
                num_cells = max(unique(regions));
                tis.cells = containers.Map('KeyType','int32','ValueType','any');
                
                % Get vertices that are not on border
                [vert_coords,vx2Cell] = tis.validate_vertices(vert_coords);
                num_vertices = size(vert_coords,1);
                vertices(1:num_vertices) = Vertex; % preallocate empties
                for i = 1:num_vertices
                    vertices(i) = Vertex(vert_coords(i,1),vert_coords(i,2));
                end
                
                % Merge vertices that re too close to each other
                tis.merge_threshold_in_px = 3;
                [vertices,vert_coords] = tis.merge_vertices(vert_coords,vertices,...
                    tis.merge_threshold_in_px);
                tis.vert_coords = vert_coords;
                tis.vertices = vertices;
                
                % Get cell-ownership of vertices via 8-connected neighbors of
                % vertices and REGIONS map
                [tis.vert_coords,vx2Cell] = tis.validate_vertices;
                
                % Instatiate valid cells
                for i = 1:num_cells
                    tis.cells(int32(i)) = ...
                        CellModel(int32(i), tis,...
                        vertices( cellfun(@(x) any(x == i),vx2Cell) ), ...
                        centroids(i,:) );
                end
                
            end
        end % Constructor
        
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
        
        % ----- Cell methods ------
        
        function tis = activateCell( tis, cellIDs )
            % Set specified cells (IDs) to "active = 1"
            % Usage: tissue = activateCell(tissue, [1 2 3]);
            %        tissue = activateCell(tissue); Activates all cells
            
            if nargin < 2,
                cellIDs = tis.cells.keys;
            end
            for i = 1:numel(cellIDs)
                tis.cells( cellIDs{i} ) = ...
                    tis.cells( cellIDs{i} ).activateCell;
            end
        end % activateCell
        
        function cells = getActiveCells(tis)
            % Return all active cells in the tissue
            % 
            % USAGE: actives = tissue.getActiveCells;
            cells = tis.cells.values;
            cells = [cells{:}];
            cells = cells( [cells.isActive] > 0 );
        end
        
        % ----- Vertex methods ------
        
        function [vert_coords,vx2Cell] = validate_vertices(tis, vert_coords)
            
            if nargin < 2, vert_coords = tis.vert_coords; end
            
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
                conn_pixels = tis.regions(x-1:x+1, y-1:y+1);
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
        
        % 		function tis = vertexCellMatrix(tis)
        %
        % 			num_cells = numel(tis.cells);
        % 			num_verts = numel(tis.vertices);
        %
        % 		end
        
        %------ Visualization -----
        function I = draw(tis)
            I = tis.regions;
        end
        
        
        
    end % Methods
    
    methods (Static)
        
        
    end
    
end
