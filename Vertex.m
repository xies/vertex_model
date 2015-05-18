classdef Vertex
    % ---- Properties ----
    %   ID
    %   x
    %   y
    %   cellIDs
    %   bondIDs
    % --- Methods ----
    %   distance - distance b/w two vertices
    %   sort - sort counter-clockwise around a given point
    %   move - move to new position
    %   rotate - rotate counter-clockwise by a given angle
    %   eq/ne - comparison done by x+y
    %   ismember - comparison by ID
    % ---Visualize ---
    %   draw
    %   line
    % 
    % xies@mit.edu March 2015
    
    properties (SetAccess = private)
        ID %key to the Map owned by Tissue.m (CellModels hold keys)
    end
    
    properties
        x
        y
        cellIDs
        bondIDs
    end
    
    methods
        
        function obj = Vertex(ID,x,y,cellIDs)
            % Vertex constructor, needs x and y coordinates as input
            if nargin > 0
                obj.ID = ID;
                obj.x = x;
                obj.y = y;
                obj.cellIDs = cellIDs;
            end
        end % Constructor
        
        % ------- Measurements --------
        function d = distance(vt1,vt2)
            % Euclidean distance
            d = sqrt(([vt1.x]-[vt2.x])^2 + ([vt1.y]-[vt2.y])^2);
        end
        
        function v_array = sortByDistance(v_array,point)
            % Sort an array based on distance wrt point
            vx = [v_array.x]; vy = [v_array.y];
            D2 = (vx - point(2)).^2 + (vy - point(1)).^2;
            [~,I] = sort(D2);
            v_array = v_array(I);
        end
        
        function v_array = sort(v_array,centroid)
            % Sort an array based on clock-wise angle wrt CENTROID
            % Should begin in the -y axis (min angle is -pi)
            vx = [v_array.x]; vy = [v_array.y];
            angles = atan2( vy-centroid(2), vx-centroid(1));
            [~,I] = sort(angles);
            v_array = v_array(I);
        end % sort
        
        % ------- Change vertex --------
        function vx = move(vx,new_pos)
            %Specify new position
            vx.x = new_pos(1); vx.y = new_pos(2);
        end
        function vt = rotate(vt,origin,theta)
            
            origin = ensure_column(origin);
            
            % Rotate a vertex with respect to an origin and by theta
            x = vt.x; y = vt.y;
            R = [cos(theta),-sin(theta);sin(theta),cos(theta)]; %2D rotation matrix
            x = x - origin(1);
            y = y - origin(2);
            
            new_vt = origin + R*[x; y];
            vt.x = new_vt(1); vt.y = new_vt(2);
            
        end
        
        % ------- Comparators --------
        function boolArr = eq(vertex_a,vertices)
            if numel(vertex_a) ~= 1 && numel(vertices) ~= 1
                error('Can''t comapre two vectors of vertices.');
            end
            boolArr = [vertex_a.x] == [vertices.x] & [vertex_a.y] == [vertices.y];
        end
        
        function boolArr = ne(vertex_a,vertices)
            if numel(vertex_a) ~= 1 && numel(vertices) ~= 1
                error('Can''t comapre two vectors of vertices.');
            end
            boolArr = [vertex_a.x] ~= [vertices.x] | [vertex_a.y] ~= [vertices.y];
        end
        
        function boolArr = ismember( tobeCompared, list)
            if isempty(list) || isempty(tobeCompared)
                boolArr = []; return;
            end
            boolArr = ismember([tobeCompared.ID], [list.ID]);
        end
        
        % ---- Visualize ----
        function draw( vArray )
            % Draw locations of vertices using SCATTER.
            x = [vArray.x]; y = [vArray.y];
            C = spring( numel(x) );
            scatter( y, x, 100, C,'filled');
        end
        function line( vArray, color)
            % Draw lines between vertices using LINE.
            if nargin < 2
                color = [0 0 1];
            end
            x = [vArray.x]; y = [vArray.y];
            line( y, x, 'Color', color);
        end
        
    end
    
end