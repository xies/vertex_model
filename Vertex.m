classdef Vertex
    % ---- Vertex ----
    % 
    % --- Properties ----
    %
    % xies@mit.edu March 2015
    
    properties
        x
        y
    end
    
    methods
        
        function obj = Vertex(x,y)
            % Vertex constructor, needs x and y coordinates as input
            if nargin > 0
                obj.x = x;
                obj.y = y;
            end
        end % Constructor
        
        % ------- Measurements --------
        function d = distance(vertex_a,vertex_b)
            % Euclidean distance
            d = sqrt((vertex_a.x-vertex_b.x)^2 + (vertex_a.y-vertex_b.y)^2);
        end
        function v_array = sort(v_array,centroid)
            % Sort an array based on clock-wise angle wrt CENTROID
            % Should begin in the -y axis (min angle is -pi)
            
            vx = [v_array.x]; vy = [v_array.y];
            angles = atan2( vy-centroid(2), vx-centroid(1));
            [~,I] = sort(angles);
            v_array = v_array(I);
            
        end
        
        % ------- Change vertex --------
        function vx = move(vx,new_pos)
            %Specify new position
            vx.x = new_pos(1); vx.y = new_pos(2);
        end
        function v = displace_by(v,vector)
            % Move vertex by a displacement vector [dx, dy]
            v.x = v.x + vector(1);
            v.y = v.y + vector(2);
        end
        
        function vx = merge(v_array)
            % Merge an array of vertices by finding mean value
            vx = Vertex( ...
                round( mean([v_array.x]) ), ...
                round( mean([v_array.y]) ) );
        end
        
        % ------- Comparator --------
        function boolArr = eq(vertex_a,vertices)
            if numel(vertex_a) ~= 1 && numel(vertices) ~= 1
                error('Can''t comapre two vectors of vertices.');
            end
            boolArr = [vertex_a.x] == [vertices.x] & [vertex_a.y] == [vertices.y];
        end
        
    end
    
end