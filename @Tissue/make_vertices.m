
function make_vertices(tis, regions)
% Checks if vertices are not beyond image border, touching at
% least one cell, and return which cells a vertex is touching
% Use only from Constructor!

% Check that tis.vertices is empty -- or else this is not
% called from the Constructor.
if tis.vertices.length > 0
    error('Don''t call this function outside of construtor');
end

vert_coords = tis.vert_coords;
num_vertices = size( vert_coords,1 );
for i = num_vertices:-1:1
    
    x = vert_coords(i,1); y = vert_coords(i,2);
    
    % Check if vx is beyond boundary, if so, chuck it and
    % move on
    if x < 1 || x > tis.Xs || y < 1 || y > tis.Ys
        vert_coords(i,:) = [];
        continue;
    end
    
    % Extract 8 connected neighbors of this vertex
    conn_pixels = regions(round(x)-1:round(x)+1, ...
        round(y)-1:round(y)+1);
    neighbor_cells = unique(conn_pixels(conn_pixels > 0));
    
    % Check that vertex has at least 1 connected cell, if no,
    % then chuck it
    if isempty(neighbor_cells)
        vert_coords(i,:) = [];
    else
        % If is not empty- then instantiate a Vertex
        % and add it to the Map (handle object so should
        % be modified in all spaces)
        tis.vertices( int32(i) ) = ...
            Vertex( int32(i), vert_coords(i,1), vert_coords(i,2), ...
            neighbor_cells);
    end
end

tis.vert_coords = vert_coords;

end % make_vertices