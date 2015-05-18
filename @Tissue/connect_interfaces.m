function connect_interfaces(tis, opt)
%CONNECT_INTERFACES
% Sets the vertex-vertex connectivity matrix according to the
% specified model configurations, and then makes the
% appropriate Interface objets.
%
% Alternatively, if 'update' is specified, will update tissue
% with correct connectivity map.
%
% USAGE:
%  tis.connect_interfaces( opt )
%  tis.connect_interfaces( opt )
%
% INPUT: tis - the tissue to be connected
%        opt - 'purse string' / 'apical' / 'both'
%              'update'
%
% @todo: implement 'apical' and 'both'

vIDsList = tis.vertices.keys; vIDsList = [vIDsList{:}];
num_vertices = numel(vIDsList);

switch opt
    % Connect the 'junctions' of cells only
    case 'purse string'
        
        % Start with brand-new hashmap
        tis.interfaces = ...
            containers.Map('KeyType','int32','ValueType','any');
        conn = zeros(num_vertices);
        ID = 1;
        for i = 1:num_vertices
            
            this_vertex = tis.getVertices( vIDsList(i) );
            
            % Get cells according to this Vertex's list
            candidateCells = tis.getCells(this_vertex.cellIDs);
            for this_cell = candidateCells
                % Find which vertices are connected through
                % this cell (should be 2 for interior cells)
                connectedVertIDs = ...
                    this_cell.getConnectedVertices( this_vertex );
                if numel(connectedVertIDs) > 2
                    keyboard
                end
                
                % Since there should be <=2 elements, for-loop
                % is fast
                for j = 1:numel(connectedVertIDs)
                    
                    that_vertex = tis.getVertices( connectedVertIDs(j) );
                    % Check for double counting: if interface
                    % already connected, then just add
                    % candicate cells to the count
                    J = find( vIDsList == connectedVertIDs(j) );
                    if conn(i,J) == 1 || conn(J,i) == 1
                        % Add cellID to edge
                        edge = tis.getInterfaceByvIDs( ...
                            [this_vertex.ID that_vertex.ID]);
                        edge.cIDs = unique([edge.cIDs this_cell.cellID]);
                        tis.interfaces(edge.ID) = edge;
                        % Add bondID to cell
                        this_cell.bondIDs = unique([this_cell.bondIDs edge.ID]);
                        tis.cells(this_cell.cellID) = this_cell;
                        continue
                    elseif conn(i,J) > 2 || conn(J,i) == 2
                        keyboard;
                    end
                    
                    % Instantiate an Interface
                    tis.interfaces(ID) = ...
                        Interface( ID,[this_vertex that_vertex], ...
                        this_cell.cellID, tis );
                    
                    % Tell vertices about interface
                    this_vertex.bondIDs = unique([this_vertex.bondIDs ID]);
                    tis.vertices( this_vertex.ID ) = this_vertex;
                    that_vertex.bondIDs = unique([that_vertex.bondIDs ID]);
                    tis.vertices( that_vertex.ID ) = that_vertex;
                    if numel( that_vertex.bondIDs) > 3
                        keyboard;
                    end
                    
                    % Tell cells about interface
                    this_cell.bondIDs = unique([this_cell.bondIDs ID]);
                    tis.cells( this_cell.cellID ) = this_cell;
                    
                    % Increment valid bIDs
                    ID = ID + 1;
                    % Keep track of which ones we've seen
                    conn(i,J) = 1; conn(J,i) = 1;
                    
                end
            end
        end
        
    case 'update'
        % Update the connectivity matrix based on the current
        % set of Interfaces.
        
        edges = tis.getInterfaces;
        %                     v = tis.getVertices;
        vIDsList = tis.vertices.keys; vIDsList = [vIDsList{:}];
        
        conn = zeros(numel(vIDsList));
        for e = edges
            I = find(vIDsList == e.vIDs(1));
            J = find(vIDsList == e.vIDs(2));
            conn( I,J ) = 1; conn( J,I ) = 1;
        end
        %                     pairs = cat(1,edges.vIDs);
        %                     conn = accumarray( pairs,1, ...
        %                         [max(vIDsList) max(vIDsList)]);
        %                     conn = conn + tril(conn.',-1);
        
    otherwise
        error('Unrecognized vertex connection option.')
end

tis.connectivity = conn;

end % connect_interfaces