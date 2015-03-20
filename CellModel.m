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
                
                if ~obj.isValid
                    error('Cannot contruct cellModel, check inputs');
                end
            end
            
        end % Constructor
        
        function flag = isValid(cellm)
            % Constructor checker -- see that everything is the right class/
            % dimensions
            flag = all(size(cellm.centroid) == [1,2]);
            % Make sure centroid given by parentTissue is also centroid of
            % the vertices
            x = [cellm.vertices.x]; y = [cellm.vertices.y];
            xplus1 = x([2:end 1]); yplus1 = y([2:end 1]);
            cx = sum( ( x+xplus1 ).*( x.*yplus1 - xplus1.*y ) )/numel(x)/cellm.area;
            cy = sum( ( y+yplus1 ).*( x.*yplus1 - xplus1.*y ) )/numel(y)/cellm.area;
            flag = flag && abs(cx - cellm.centroid(1)) <= 1 && abs(cy - cellm.centroid(2)) <= 1;
            if flag < 1, keyboard; end
        end
        
        % --------- Measurements ---------
        
        function a = get_anisotropy(cellm)
            I = cellm.draw_smallMask;
            a = regionprops(I,'MajorAxisLength','MinorAxisLength');
            a = a.MinorAxisLength/a.MajorAxisLength;
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
        
        % --------- Vertex methods --------------
        
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
            
        end
        
        % --------- Set contractility --------------
        function cell = activateCell(cell)
            cell.isActive = 1;
        end
        
        % ---------  Cell-cell connectivity --------
        
        function cells = get_neighbors(this_cell,n)
            % GET_NEIGHBORS
            %
            % neighbors = this_cell.get_neighbors(n);
            
            if nargin < 2, n = 1; end % By default give first-order neighbors
            tissue = this_cell.parent;
            cells = tissue.get_neighbors(this_cell.cellID,n);
        end % get_neighbors
        
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
            Xs = max(vx) + 1; Ys = max(vy) + 1;
            image = poly2mask(vy,vx,Xs,Ys);
        end
        
    end % End methods
    
    methods (Static)
        
    end
    
end
