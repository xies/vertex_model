classdef CellModel
    % ---- CellModel ----
    % Base class for vertex model; contains one cell in a tissue at one
    % single time point
    %
    % --- Properties ----
    %   centroid
    %   vertices
    %   area
    %   contractility
    %   anisotropy
    %
    %   cellID
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
        vertices	% vertices, y-coordinates
        
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
            % Parent from Tissue.m
            
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
                
%                 if ~obj.isValid
%                     error('Cannot contruct cellModel, check inputs');
%                 end
            end
            
        end % Constructor
        
        function flag = isValid(cellm)
            % Constructor checker -- see that everything is the right class/
            % dimensions
            flag = all(size(cellm.centroid) == [1,2]);
            % Make sure centroid given by parentTissue is also centroid of
            % the vertices
            ct = cellm.get_centroid;
            flag = flag && abs(ct(1) - cellm.centroid(1)) <= 1;
            flag = flag && abs(ct(2) - cellm.centroid(2)) <= 1;
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
            % @todo: Fix this situation!!!
            I = cellm.draw_smallMask;
            a = regionprops(I,'MajorAxisLength','MinorAxisLength');
            if numel(a) == 1
                a = a.MinorAxisLength/a.MajorAxisLength;
            else
                a = NaN;
            end
        end
        function a = get_area(cellm)
            x = [cellm.vertices.x]; y = [cellm.vertices.y];
            xplus1 = x([2:end 1]); yplus1 = y([2:end 1]);
            a = sum( (x.*yplus1) - (y.*xplus1) )/2;
        end
        function p = get_perimeter(cellm)
            x = [cellm.vertices.x]; y = [cellm.vertices.y];
            xplus1 = x([2:end 1]); yplus1 = y([2:end 1]);
            p = sum( sqrt((x-xplus1).^2 + (y-yplus1).^2) );
        end
        function centroid = get_centroid(cellm)
            % Calculate centroid from vertices (x,y)
            x = [cellm.vertices.x]; y = [cellm.vertices.y];
            xplus1 = x([2:end 1]); yplus1 = y([2:end 1]);
            cx = sum( ( x+xplus1 ).*( x.*yplus1 - xplus1.*y ) )/numel(x)/cellm.area;
            cy = sum( ( y+yplus1 ).*( x.*yplus1 - xplus1.*y ) )/numel(y)/cellm.area;
            centroid = [cx, cy];
        end
        
        % ------- Cell set methods -------
        
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
        
        % --------- Vertex set methods --------------
        
        function c_array = sort(c_array,centroid)
            % Sort an array based on clock-wise angle wrt CENTROID
            % Should begin in the -y axis (min angle is -pi)
            
            ct = cat(1,c_array.centroid);
            angles = atan2( ct(:,1)-centroid(2), ct(:,2)-centroid(1));
            [~,I] = sort(angles);
            c_array = c_array(I);
            
        end
        
        function c = moveVertex(c,vert,new_vcoord)
            % Move vertex contained in C
            % USAGE: c = c.moveVertex(vertex,new_vcoord);
            
            if ~any(c.vertices == vert)
                error('vertex does not belong to this cell');
            end
            
            I = find(c.vertices == vert);
            c.vertices(I) = c.vertices(I).move( new_vcoord );
            c = c.updateCell;
            
        end
        
        % --------- Set contractility --------------
        function cell = activateCell(cell)
            cell.isActive = 1;
        end
        function cell = deactivateCell(cell)
            cell.isActive = 0;
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
