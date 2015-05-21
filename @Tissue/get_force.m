function F = get_force(tis)
% GET_FORCE
%
% Returns the forces as given by -grad(E) / drag_coeff of
% current tissue configuration.
%
% USAGE: V = tis.get_velocities;
%
% Currently implements, for all cells c, vertices i (and its neighbor j):
%
%   - area elasticity:
%       Fa = -sum_c( k*(Ac-A0)*grad_i(Ac) )
%   - perimeter elasticity
%       Fp = -sum_c( -k*(Pc-P0)*grad_i(Pc)
%   - line tension:
%       Fl = -sum_i( sum_j(sij*grad_i(Dij)) )
%
% The gradient terms are implemented as Nagai (2000).

% Grab data
vcoords = tis.vert_coords;
D = tis.interVertDist;

% Grap tissue-constant parameters
length_scale = tis.parameters.lengthScale;

% Initialize
num_verts = size(vcoords,1);
F = zeros( size(vcoords) );

vIDList = tis.vertices.keys;
vIDList = [vIDList{:}];
% For now, loop through; vectorize later?
for i = 1:num_verts
    
    if tis.parameters.fixed_verts(i)
        continue;
    end
    
    % Find connected vertices
    vi = tis.vertices(vIDList(i));
    num_neighbors = numel( vi.cellIDs );
    
    % Only give nonzero velocities for vertices w/ more than 2
    % cell neighbors (non-fixed cells)
    if num_neighbors > 2
        
        % Find conn cells (this is much faster)
        neighbors(1:num_neighbors) = CellModel(); % initialize
        edges(1:num_neighbors) = Interface(); % initialize
        for j = 1:num_neighbors
            edges(j) = tis.interfaces( vi.bondIDs(j) );
            neighbors(j) = tis.cells(vi.cellIDs(j));
        end
        
        % Intialize
        line_tension_term = [0 0];
        area_elastic_term = [0 0];
        perim_elastic_term = [0 0];
        active_contraction_term = [0 0];
        
        % Generate line tension term in direction of the
        % inter-vertex vector
        % Since loop is only a few elements long, it's faster
        % than vectorized version
        for j = 1:numel(edges)
            
            other_vert = edges(j).vIDs;
            other_vert = other_vert(other_vert ~= vi.ID);
            other_vert = vIDList == other_vert;
            line_tension_term = line_tension_term - ...
                edges(j).tension ...
                * (vcoords(i,:) - vcoords(other_vert,:)) ...
                / D(i,other_vert);
        end
        
        % Go through all CELLS associated with current vertex,
        % and calculate cell elasticity.
        % NOTE that we can't just go through "edges" themselves
        for this_cell = neighbors
            
            % Need to sort vertices counter-clockwise
            sortedVt = tis.getVertices( this_cell.vIDs );
            sortedVt = sortedVt.sortClockwise( this_cell.centroid );
            
            % Find this current vertex and its immediate
            % neighboring vertices
            I = find( vi == sortedVt );
            I = wrap( [I-1 I I + 1] , numel(sortedVt)); % circularly index
            r = [sortedVt(I).x ; sortedVt(I).y];
            
            R = [0 -1; 1 0]; % pi/2 rotation matrix
            v = ( R*(r(:,1) - r(:,3)) )'; % grad(A)
            
            ua = r(:,1) - r(:,2); ub = r(:,3) - r(:,2);
            u = ua / norm(ub) + ub /norm(ub); % Direction of grad(P)
            
            % Area elasticity
            area_elastic_term = area_elastic_term ...
                - 2 * (this_cell.area - this_cell.targetArea) ...
                * v * this_cell.areaElasticity / length_scale ^2;
            
            % Perimeter elasticity
            perim_elastic_term = perim_elastic_term ...
                + 2 * (this_cell.perimeter) * u' ...
                * this_cell.perimElasticity / length_scale;
            
            % Active contraction -- same form as area
            % elasticity but sets A0 = 0 and uses different
            % coefficient
            active_contraction_term = active_contraction_term ...
                - 2 * this_cell.area * this_cell.contractility ...
                * v / length_scale^2;
            
        end
        
        F(i,:) = line_tension_term + area_elastic_term + ...
            perim_elastic_term + active_contraction_term;
        
        %                     tis.draw('showVectors',{V(i,:),i},'showActive');
        %                     drawnow
        
    end
    
end

end % get_forces