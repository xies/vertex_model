classdef CellModel
    % ---- CellModel ----
    % Base class for vertex model; contains one cell in a tissue at one
    % single time point
    %
    % --- Properties ----
    %
    %   (private)
    %      cellID - key used in parent Map (int32)
    %      parentTissue - Tissue object
    %   (public)
    %      centroid - centroid coord
    %      vertices - array of Vertex objects
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
        parentTissue    % the Tissue parent
        
    end % Private properties (can't be set by outside methods)
    
    properties
        
        % measurements
        centroid	% centroid of cell
        vertices	% vertices, xy-coordinates
        
        contractility % current contractility
        anisotropy  % anisotropy (shape)
        area		% area
        perimeter	% perimeter
        % simulation properties
        isActive    % active for contractility
        
    end % Public properties
    
    methods
        
        function obj = CellModel(ID,parent,vt,ct)
            % CellModel - Contructs object of CellModel class.
            % Parent is Tissue.m
            % 
            % USAGE: cellm = CellModel( ID, parent, verts, centroid)
            %
            % INPUT: ID - int32 cellID
            %        parent - parent Tissue.m
            %        verts - array of Vertex objects
            %        centroid - (x,y) coord of centroid
            
            if nargin > 0
                
                % Inherit properties from tissue
                obj.cellID = ID;
                obj.centroid = ct;
                obj.vertices = sort(vt,ct); %sort vertices so they're clockwise (implemented in Vertex.m)
                obj.parentTissue = parent;
                
                % Make measurements based on current setting
                obj.anisotropy = obj.get_anisotropy;
                obj.perimeter = obj.get_perimeter;
                obj.area = obj.get_area;
                
                % Put in default values
                obj.contractility = 0;
                obj.isActive = 0;
                
                if ~obj.isValid
                    error('Cannot contruct cellModel, check inputs');
                end
            end
            
        end % Constructor
        
        function flag = isValid(cellm)
            % Constructor checker -- see that everything is the right class/
            % dimensions and current centroid is the centroid of current
            % set of vertices
            %
            % USAGE: flag = isValid( cellm );
            
            flag = all(size(cellm.centroid) == [1,2]);
            % Make sure centroid given by parentTissue is also centroid of
            % the vertices
            ct = cellm.get_centroid;
            flag = flag && abs(ct(1) - cellm.centroid(1)) <= 1;
            flag = flag && abs(ct(2) - cellm.centroid(2)) <= 1;
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
        
        function a = get_anisotropy(cellm)
            % Measures anisotropy of current cell by drawing a mask of it
            % and running regionprops
            % 
            % USAGE: a = get_anisotropy(cellm)
            %
            % @todo: Fix this situation where object breaks by not using
            % the binary mask but vertices themselves
            
            I = cellm.draw_smallMask;
            a = 0;
%             a = regionprops(I,'MajorAxisLength','MinorAxisLength');
%             if numel(a) == 1
%                 a = a.MinorAxisLength/a.MajorAxisLength;
%             else
%                 a = NaN;
%             end
        end
        function a = get_area(cellm)
            % Calculates area from the vertex positions directly
            % USAGE: a = get_area(cellm)
            x = [cellm.vertices.x]; y = [cellm.vertices.y];
            xplus1 = x([2:end 1]); yplus1 = y([2:end 1]);
            a = sum( (x.*yplus1) - (y.*xplus1) )/2;
        end
        function p = get_perimeter(cellm)
            % Calculates perimeter from the vertex positions directly
            % USAGE: p = get_perimeter(cellm)
            x = [cellm.vertices.x]; y = [cellm.vertices.y];
            xplus1 = x([2:end 1]); yplus1 = y([2:end 1]);
            p = sum( sqrt((x-xplus1).^2 + (y-yplus1).^2) );
        end
        function centroid = get_centroid(cellm)
            % Calculates centroid from the vertex positions directly
            % USAGE: ct = get_centroid(cellm)
            x = [cellm.vertices.x]; y = [cellm.vertices.y];
            G = polygeom(x,y);
%             xplus1 = x([2:end 1]); yplus1 = y([2:end 1]);
%             cx = sum( ( x+xplus1 ).*( x.*yplus1 - xplus1.*y ) )/numel(x)/cellm.area;
%             cy = sum( ( y+yplus1 ).*( x.*yplus1 - xplus1.*y ) )/numel(y)/cellm.area;
            centroid = [G(2), G(3)];
        end
        
        % ------- Cell set methods -------
        
        function cell = activateCell(cell), cell.isActive = 1; end
        function cell = deactivateCell(cell), cell.isActive = 0; end
        
        function cellm = setContractility( cellm, C)
            % Sets the active contractility coefficient of a cell
            cellm.contractility = C;
        end
        
        function cellm = updateCell(cellm)
            % Updates the centroid to the current centroid, area,
            % perimeter, anisotropy based on the new list of vertices.
            % 
            % USAGE: cell = updateCentroid(cell)
            ct = cellm.get_centroid;
            cellm.centroid = ct;
            cellm.area = cellm.get_area;
            cellm.perimeter = cellm.get_perimeter;
            cellm.anisotropy = cellm.get_anisotropy;
            
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
        
        
        % --------- Vertex set methods --------------
        
        function c = moveVertex(c,vert,new_vcoord)
            % Move vertex contained in input cell and the update all
            % properties of the cell to match new vertices
            %
            % USAGE: c = c.moveVertex(vertex,new_vcoord);
            
            if ~any(c.vertices == vert)
                error('vertex does not belong to this cell');
            end
            
            I = find(c.vertices == vert);
            c.vertices(I) = c.vertices(I).move( new_vcoord );
            c = c.updateCell;
            
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
            % USAGE: neighborVs = getConnectedVertices(cellm,v)
            
            I = find(cellm.vertices == v);
            if isempty(I), neighborVerts = []; return; end
            
            if I == 1
                neighborVerts = [cellm.vertices(end), cellm.vertices(2)];
            elseif I == numel( cellm.vertices )
                neighborVerts = [cellm.vertices(end-1), cellm.vertices(1)];
            else
                neighborVerts = [cellm.vertices(I-1), cellm.vertices(I+1)];
            end
                
        end
        
        % --------- Visualize ---------
        
        function image = draw(cellm)
            % Draws a binary outline of the cell
            vx = [cellm.vertices.x]; vy = [cellm.vertices.y];
            [xe,ye] = poly2edge(vx,vy);
            Xs = cellm.parentTissue.Xs;
            Ys = cellm.parentTissue.Ys;
            image = accumarray(round([xe,ye]), 1, [Xs Ys], @max);
        end
        
        function image = drawMask(cellm)
            % Use POLY2MASK to draw a mask of the current cell
            v = cellm.vertices;
            Xs = cellm.parentTissue.Xs;
            Ys = cellm.parentTissue.Ys;
            image = poly2mask([v.y],[v.x],Xs,Ys);
        end
        function image = draw_smallMask(cellm)
            % Use POLY2MASK to draw a mask of the current cell--only around
            % bounding box of cell
            v = cellm.vertices;
            vx = [v.x]; vy = [v.y];
            vx = vx - min(vx) + 1; vy = vy - min(vy) + 1;
            Xs = round(max(vx)) + 1; Ys = round(max(vy)) + 1;
            image = poly2mask(vy,vx,Xs,Ys);
        end
        
    end % End methods
    
end
