classdef CellModel
    % ---- CellModel ----
    % Base class for vertex model; contains one cell in a tissue at one
    % single time point
    %
    % --- Properties ----
    %
    %   (private)
    %      cellID - key used in parent Map (int32)
    %   (public)
    %      centroid - centroid coord
    %      vIDs - array of indices to Vertex objects
    %      bondIDs - array of indices to Interface objects
    %      area     - cell area
    %      perimeter  - perimeter
    %      anisotropy - ellipsoid anisotropy
    %      energy - sum of contractile, area and perim elasticity
    %      --- Parameters ---
    %      targetArea - A0
    %      areaElasticity
    %      perimElasticity
    %      
    %      contractility - current active force amt
    %      isActive - active/contracting flag
    %
    % --- Methods ----
    %      CellModel - constructor
    %      isValid - consistency checker
    % 
    %   --- Measurements ---
    %       get_area - get area from vertex positions
    %       get_perimeter - get perimeter from vertex positions
    %       get_centroid - get centroid from vertex positions
    %       get_anisotropy - get area from vertex-derived mask  (@todo work with non-binary)
    %       get_distance_to - get distance of cell centroid to another
    %           given point
    %       contains - whether a point is contained in a cell
    %       get_vertices - returns the set of Vertex models touching this
    %          cell
    %       get_bonds - returns the set of Interface models touching this
    %          cell
    %
    %   --- Cell set methods ---
    %       activateCell - set cell to active
    %       deactivateCell -  set cell to not active
    %       setContractility - sets the contractility of a cell
    %       updateCell - update all properties of CellModel to match
    %          current set of vertices
    %       sort - sorts an array of cells in counter-clockwise w.r.t. an
    %          input point
    %       
    %   --- Vertex set methods ---
    %       moveVertex - moves a cell's vertices to new positions and
    %          update cell (uses updateCell)
    %
    %   --- Connectivity ---
    %       connected - sees if two vertices are edge-connected within an input
    %          cell
    %       getConnectedVertices - returns all edge-connected vertices
    %          (should be 2 of them for all inputs) of an input vertex as
    %          seen by a single cell
    %
    %   --- Visualize ---
    %       draw - draws the outline of the cell in binary image
    %       drawMask - draws the mask of cell in binary image
    %       draw_smallMask - draws a small mask of cell
    %
    % xies@mit.edu March 2015
    
    properties (SetAccess= private)
        
        % IDs
        cellID	% cellID in EDGE
        
    end % Private properties (can't be set by outside methods)
    
    properties
        
        % measurements
        centroid	% centroid of cell
        vIDs	% IDs to Vertex objects
        bondIDs    % IDs to interface objects
        
        contractility % current contractility
        anisotropy  % anisotropy (shape)
        area		% area
        perimeter	% perimeter
        energy
        
        % simulation properties
        isActive    % active for contractility
        areaElasticity
        perimElasticity
        targetArea
        
    end % Public properties
    
    methods
        
        function obj = CellModel(ID,parent,vIDs,ct)
            % CellModel - Contructs object of CellModel class.
            % Parent is Tissue.m
            % 
            % USAGE: cellm = CellModel( ID, parent, verts, centroid)
            %
            % INPUT: ID - int32 cellID
            %        parent - parent Tissue.m
            %        verts - array of keys to tis.vertices map
            %        centroid - (x,y) coord of centroid
            
            if nargin > 0
                
                % Inherit properties from tissue
                obj.cellID = ID;
                obj.centroid = ct;
                
                % Sort all vertices
                vts = parent.getVertices( vIDs );
                vts = vts.sortClockwise( obj.centroid );
                obj.vIDs = [vts.ID];

                % Put in default values
                obj.contractility = 0;
                obj.isActive = 0;
                
                % -- set measurements ---
                obj.area = obj.get_area(parent);
                obj.anisotropy = obj.get_anisotropy(parent);
                obj.perimeter = obj.get_perimeter(parent);
                
                if ~obj.isValid(parent)
                    error('Cannot contruct cellModel, check inputs');
                end
            end
            
        end % Constructor
        
        function flag = isValid(cellm,tis)
            % Constructor checker -- see that everything is the right class/
            % dimensions and current centroid is the centroid of current
            % set of vertices
            %
            % USAGE: flag = isValid( cellm );
            
            flag = all(size(cellm.centroid) == [1,2]);
            % Make sure centroid given by parentTissue is also centroid of
            % the vertices
            ct = cellm.get_centroid(tis);
            flag = flag && abs(ct(1) - cellm.centroid(1)) <= 2;
            flag = flag && abs(ct(2) - cellm.centroid(2)) <= 2;
            if ~flag, keyboard; end
        end % isValid
        
        function flags = eq( this_cell, cellArray)
            % Evaluates whether a cell is equal to any element of an array
            % of cells
            
            ct = cat(1,cellArray.centroid);
            flags = any( this_cell.centroid(1) == ct(:,1) ); % cx
            flags = flags & any( this_cell.centroid(2) == ct(:,2) ); % cy
            
        end % eq
        
        % --------- Measurements ---------
        
        function E = get_energy(cellm,tis)
            %GET_ENERGY
            % Caculates the current reduced energy associated with this cell
            % E = area_elastic + perim_minization + contractility
            
            if tis.parameters.dimensionless
                sigma = tis.parameters.forceScale;
                lambda = tis.parameters.lengthScale;
            else
                sigma = 1;
                lambda = 1;
            end
                E = cellm.areaElasticity*(cellm.area - cellm.targetArea)^2 / 2;
                E = E + cellm.perimElasticity*(cellm.perimeter)^2 / 2;
                E = E + cellm.contractility * cellm.area ^ 2;
                
                E = E/sigma/lambda;
        end % get_energy
        
        function a = get_anisotropy(cellm,tis)
            %GET_ANISOTROPY
            % Measures anisotropy of current cell by drawing a mask of it
            % and running regionprops
            % 
            % USAGE: a = get_anisotropy(cellm,tis)
            %
            % @todo: Fix this situation where object breaks by not using
            % the binary mask but vertices themselves
            a = 1;
%             vt = cellm.get_vertices(tis);
%             f = fit_ellipse([vt.y],[vt.x]);
%             if ~isempty(f)
%                 a = f.short_axis / f.long_axis;
%             else
%                 a = NaN;
%             end
        end % get_anisotropy
        
        function a = get_area(cellm,tis)
            %GET_AREA
            % Calculates area from the vertex positions directly
            %
            % USAGE: a = get_area(cellm,tis)
            %
            % See also: POLYAREA
            
            vt = tis.getVertices( cellm.vIDs );
%             vt = sort(vt, cellm.centroid); % sort counter-clockwise
            x = [vt.x]; y = [vt.y];
            if ~isempty(tis.parameters)
                l = tis.parameters.um_per_px;
            else
                l = 1;
            end
            a = polyarea( x, y) * l^2;
            
        end % get_area
        
        function p = get_perimeter(cellm,tis)
            %GET_PERIMETER
            % Calculates perimeter from the vertex positions directly
            %
            % USAGE: p = get_perimeter(cellm,tis)
            vt = tis.getVertices( cellm.vIDs );
%             vt = vt.sort(cellm.centroid);
            x = [vt.x]; y = [vt.y];
            xplus1 = x([2:end 1]); yplus1 = y([2:end 1]);
            p = sum( sqrt((x-xplus1).^2 + (y-yplus1).^2) );
        end % get_perimeter
        
        function centroid = get_centroid(cellm,tis)
            %GET_CENTROID
            % Calculates centroid from the vertex positions directly
            % USAGE: ct = get_centroid(cellm,tis)
            
            vt = tis.getVertices( cellm.vIDs );
            vt = sortClockwise(vt, cellm.centroid); % sort counter-clockwise
            x = [vt.x]; y = [vt.y];
            G = polygeom(x,y);
            centroid = [G(2), G(3)];
        end % get_centroid
        
        function D = get_distance_to(cellm,pt)
            %GET_DISTANCE_TO
            % Get distance of cell centroid to external points
            % USAGE: D = cellm.get_distance_to(points)
            ct = cellm.centroid;
            D = sqrt( (ct(1) - pt(:,1)).^2 + (ct(2) - pt (:,2)).^2 );
        end % get_distance_to
        
        function theta = get_angle(cells)
            %GET_ANGLE
            % Get the angle that the centroids of two cellmodels make to
            % the horizontal x-axis
            % 
            % USAGE: theta = cells.get_angle;
            % -pi < theta < pi
            if numel(cells) ~= 2, error('Need precisely 2 cells'); end
            ct = cat(1,cells.centroid);
            theta = atan2( diff(ct(:,1)),diff(ct(:,2)) );
        end
        
        function vts = get_vertices(cellm,tis)
            vts = tis.getVertices( cellm.vIDs );
        end
        
        function bds = get_bonds(cellm,tis)
            bds = tis.getInterfaces( cellm.bondIDs );
        end
        
        function in = contains(cellm,point,tis)
            %contains
            % Test for whether given point is inside given cell.
            %
            % USAGE: flag = cellm.contains( point, tis )
            
            vt = tis.getVertices( cellm.vIDs );
            in = inpolygon( point(2),point(1),[vt.y],[vt.x]);
            
        end
        
        % ------- Cell set methods -------
        
        function cell = activateCell(cell), cell.isActive = 1; end
        function cell = deactivateCell(cell)
            % Deactivate the active flag and sets contractility to 0
            cell.isActive = 0;
            cell.contractility = 0;
        end
        
        function cellm = setContractility( cellm, C)
            % Sets the active contractility coefficient of a cell
            cellm.contractility = C;
        end
        
        function cellm = updateCell(cellm,tis)
            % Updates the CellModel to the current centroid, area,
            % perimeter, anisotropy based on the new list of vertices.
            % 
            % USAGE: cell = cell.updateCell(tis)
            
            % Make sure cell vts are sorted according to edges!
            vts = tis.getVertices(cellm.vIDs);
            vts_old = vts;
            % Start with the most "bottom" vertex and move left + up
            [~,I] = min( [vts.x] ); vts(1) = vts(I);
            backtrace = [];
            for i = 2:numel(vts)
                % Use which edge we were just at (backtrace) to trace
                % clockwise
                vts(i) = vts(i-1).next(cellm,tis,backtrace);
                backtrace = vts(i-1);
            end
            
            if unique([vts.ID]) ~= unique([vts_old.ID])
                keyboard;
            end
            
            cellm.vIDs = [vts.ID];
            % Check that we have positive area -- if not, need to reverse
            cellm.area = cellm.get_area( tis );
            if cellm.area < 0
                cellm.vIDs = cellm.vIDs(end:-1:1);
                cellm.area = cellm.get_area( tis );
            end
            % Fill in other info
            cellm.centroid = cellm.get_centroid( tis );
            cellm.perimeter = cellm.get_perimeter( tis );
            cellm.anisotropy = cellm.get_anisotropy( tis );
            cellm.energy = cellm.get_energy(tis);
            
        end % updatedCell
        
        function c_array = sortByDistance(c_array,point)
            % Sort an array based on distance wrt a POINT
            % 
            % USAGE: cells_array = sortByDistance(cell_array, point)
            
            ct = cat(1,c_array.centroid);
            D2 = (ct(:,1)-point(2)).^2 + (ct(:,2)-point(1)).^2;
            [~,I] = sort(D2);
            c_array = c_array(I);
        end
        
        function c_array = sortByAngle(c_array,point)
            % Sort an array based on clock-wise angle wrt a POINT
            % 
            % USAGE: cells_array = sortByAngle(cell_array, point)
            
            ct = cat(1,c_array.centroid);
            angles = atan2( ct(:,1)-point(2), ct(:,2)-point(1));
            [~,I] = sort(angles);
            c_array = c_array(I);
            
        end % sort
        
        function c_array = sortByID( c_array )
            % Sort a CellModel array by ID.
            IDs = [c_array.cellID];
            [~,I] = sort(IDs);
            c_array = c_array(I);
        end
        
        % ---------  Connectivity --------
          
        function flag = connected(cells,vertA,vertB)
            % Returns if two vertices are connected by an edge in this
            % cell. Calculates whether vertA and vertB are both in the
            % cell, and also are next to each other in the sorted list.
            % 
            % USAGE: flag = cellm.connected(vertA, vertB)
            
            flag = 0;
            verts = cells.vertices;
            % Check first that both are in cell
            if all( verts ~= vertA) || all( verts ~= vertB)
                return;
            end
            % Then check for vertex order
            vertsplus1 = verts( [2:end 1] );
            for i = 1:numel(verts)
                if verts(i) == vertA && vertsplus1(i) == vertB
                    flag = 1; return;
                elseif verts(i) == vertB && vertsplus1(i) == vertA
                    flag = 1; return;
                end
            end
            
        end % connected
        
        function neighborVerts = getConnectedVertices( cellm, v)
            % Returns all the edge-connected vertex of v, as seen by cellm
            % USAGE: neighborVs = getConnectedVertices(cellm,tis,v)
            
            I = find(cellm.vIDs == v.ID);
            if isempty(I), neighborVerts = []; return; end
            
            if I == 1
                neighborVerts = [cellm.vIDs(end), cellm.vIDs(2)];
            elseif I == numel( cellm.vIDs )
                neighborVerts = [cellm.vIDs(end-1), cellm.vIDs(1)];
            else
                neighborVerts = [cellm.vIDs(I-1), cellm.vIDs(I+1)];
            end
                
        end % getConnectedVertices
        
        % --------- Visualize ---------
        
        function image = draw(cellm,tis)
            % Draws a binary outline of the cell
            v = tis.getVertices( cellm.vIDs );
%             v = v.sort(cellm.centroid);
            vx = [v.x]; vy = [v.y];
            [xe,ye] = poly2edge(vx,vy);
            Xs = tis.Xs; Ys = tis.Ys;
            image = accumarray(round([xe,ye]), 1, [Xs Ys], @max);
%             if cellm.cellID == 41, keyboard; end
        end % draw
        
        function image = drawMask(cellm,tis)
            % Use POLY2MASK to draw a mask of the current cell
            v = tis.getVertices( cellm.vIDs );
%             v = v.sort(cellm.centroid);
            Xs = tis.Xs; Ys = tis.Ys;
            image = poly2mask([v.y],[v.x],Xs,Ys);
            
        end % drawMask
        
        function image = draw_smallMask(cellm,tis)
            % Use POLY2MASK to draw a mask of the current cell--only around
            % bounding box of cell
            v = tis.getVertices( cellm.vIDs );
            vx = [v.x]; vy = [v.y];
            vx = vx - min(vx) + 1; vy = vy - min(vy) + 1;
            Xs = round(max(vx)) + 1; Ys = round(max(vy)) + 1;
            image = poly2mask(vy,vx,Xs,Ys);
            
        end % draw_smallMask
        
    end % End methods
    
end
