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
    %      area     - cell area
    %      perimeter  - perimeter
    %      anisotropy - ellipsoid anisotropy
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
        
        contractility % current contractility
        anisotropy  % anisotropy (shape)
        area		% area
        perimeter	% perimeter
        % simulation properties
        isActive    % active for contractility
        
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
                vts = vts.sort( obj.centroid );
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
        end
        
        function flags = eq( this_cell, cellArray)
            % Evaluates whether a cell is equal to any element of an array
            % of cells
            
%             num_comp = numel(cellArray);
            ct = cat(1,cellArray.centroid);
            flags = any( this_cell.centroid(1) == ct(:,1) ); % cx
            flags = flags & any( this_cell.centroid(2) == ct(:,2) ); % cy
            
%             vt = {cellArray.vertices};
%             cellfun(@(x)
            
        end
        
        % --------- Measurements ---------
        
        function a = get_anisotropy(cellm,tis)
            % Measures anisotropy of current cell by drawing a mask of it
            % and running regionprops
            % 
            % USAGE: a = get_anisotropy(cellm,tis)
            %
            % @todo: Fix this situation where object breaks by not using
            % the binary mask but vertices themselves
            a = 1;
        end % get_anisotropy
        
        function a = get_area(cellm,tis)
            % Calculates area from the vertex positions directly
            %
            % USAGE: a = get_area(cellm,tis)
            %
            % See also: POLYAREA
            
            vt = tis.getVertices( cellm.vIDs );
%             vt = sort(vt, cellm.centroid); % sort counter-clockwise
            x = [vt.x]; y = [vt.y];
            a = polyarea( x, y);
        end % get_area
        
        function p = get_perimeter(cellm,tis)
            % Calculates perimeter from the vertex positions directly
            %
            % USAGE: p = get_perimeter(cellm,tis)
            
            vt = tis.getVertices( cellm.vIDs );
            x = [vt.x]; y = [vt.y];
            xplus1 = x([2:end 1]); yplus1 = y([2:end 1]);
            p = sum( sqrt((x-xplus1).^2 + (y-yplus1).^2) );
        end % get_perimeter
        
        function centroid = get_centroid(cellm,tis)
            % Calculates centroid from the vertex positions directly
            % USAGE: ct = get_centroid(cellm,tis)
            
            vt = tis.getVertices( cellm.vIDs );
%             vt = sort(vt, cellm.centroid); % sort counter-clockwis
            x = [vt.x]; y = [vt.y];
            G = polygeom(x,y);
            centroid = [G(2), G(3)];
        end % get_centroid
        
        % ------- Cell set methods -------
        
        function cell = activateCell(cell), cell.isActive = 1; end
        function cell = deactivateCell(cell), cell.isActive = 0; end
        
        function cellm = setContractility( cellm, C)
            % Sets the active contractility coefficient of a cell
            cellm.contractility = C;
        end
        
        function cellm = updateCell(cellm,tis)
            % Updates the CellModel to the current centroid, area,
            % perimeter, anisotropy based on the new list of vertices.
            % 
            % USAGE: cell = cell.updateCell(tis)
            cellm.centroid = cellm.get_centroid( tis );
            cellm.area = cellm.get_area( tis );
            cellm.perimeter = cellm.get_perimeter( tis );
            cellm.anisotropy = cellm.get_anisotropy( tis );
            
        end
        
        function c_array = sort(c_array,point)
            % Sort an array based on clock-wise angle wrt a POINT
            % 
            % USAGE: cells_array = sort(cell_array, point)
            
            ct = cat(1,c_array.centroid);
            angles = atan2( ct(:,1)-point(2), ct(:,2)-point(1));
            [~,I] = sort(angles);
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
        end
        
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
                
        end
        
        % --------- Visualize ---------
        
        function image = draw(cellm,tis)
            % Draws a binary outline of the cell
            v = tis.getVertices( cellm.vIDs );
%             v = sort(v, cellm.centroid); % sort counter-clockwis
            vx = [v.x]; vy = [v.y];
            [xe,ye] = poly2edge(vx,vy);
            Xs = tis.Xs; Ys = tis.Ys;
            image = accumarray(round([xe,ye]), 1, [Xs Ys], @max);
        end
        
        function image = drawMask(cellm,tis)
            % Use POLY2MASK to draw a mask of the current cell
            v = tis.getVertices( cellm.vIDs );
%             v = sort(v, cellm.centroid); % sort counter-clockwis
            Xs = tis.Xs; Ys = tis.Ys;
            image = poly2mask([v.y],[v.x],Xs,Ys);
        end
        function image = draw_smallMask(cellm,tis)
            % Use POLY2MASK to draw a mask of the current cell--only around
            % bounding box of cell
            v = tis.getVertices( cellm.vIDs );
%             v = sort(v, cellm.centroid); % sort counter-clockwis
            vx = [v.x]; vy = [v.y];
            vx = vx - min(vx) + 1; vy = vy - min(vy) + 1;
            Xs = round(max(vx)) + 1; Ys = round(max(vy)) + 1;
            image = poly2mask(vy,vx,Xs,Ys);
        end
        
    end % End methods
    
end
