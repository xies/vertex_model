classdef Tissue
    % TISSUE
    %
    % ---- Properties ----
    %   vert_coords - Nv x 2 array of vertex locations
    %   vertices % hashmap of Vertex objects keyed by ID
    %   interfaces % hashmap of Interface.m objects
    %   cells - container.Map hashmap of CellModels keyed by cellID
    %   connectivity - Nv x Nv adjacency matrix of how vertices are
    %          connected
    %   interVertDist - Nv x Nv distance matrix of connected vertices
    %
    %   parameters - see setParameters for details
    %   contraction_param - contractile parameters
    %
    %   Ys - tissue image size (for display and initiation only)
    %   Xs - tissue image size
    %   merge_threshold_in_px - merge all vertices closer than this threshold
    %   t - time stamp
    %
    % ---- Methods ----
    % 
    %   Tissue - constructor
    %   isValid - consistency checker
    % 
    %   --- Energy methods ---
    %       get_energy - calculate the pot. energy of current config
    %       get_velocity - calculates the velocity on each vertex
    %
    %   --- Measurements ---
    %       getAreas - return area as function of sim time
    %       getContractility - contractility
    %       getCentroidX - x-coordinates
    %       getCentroidY - y-coordinates
    %       getTime - simulation time
    %   
    %   --- Simulation housekeeping methods ---
    %       evolve - make new configuration via updating vert_coords
    %       serParameters - sets the parameters and connects vertices
    %       activateCell - set certain cells to be active
    %       deactivateCell  - set certain cells to be not active
    %       deactivateBorder - set the outermost cells to be not active
    %
    %   --- Vertex-vertex connectivity ---
    %       connectVertices - sets the 'edges'/'interfaces' of the model
    %             upon which energy/force is calculated
    %
    %   --- Vertex-cell connectivity ---
    %       cellsContainingVertex - returns all cells that contain input
    %              vertex
    %
    %   --- Cell-cell connectivity ---
    %       connected - whether two cells share one vertex
    %       allNthOrderNeighbors - returns the filled-in "halo" of
    %              n-th order neighbors of input cell
    %       neighborsOfCell - returns just the ring of n-th order neighbors
    %              of input cell
    %       numCellNeighbors - returns # of immediate neighbors of input
    %              cell
    % 
    %   --- Cell/vertex handling ---
    %       getActiveCells - return CellModels that are active
    %       getInactiveCells - return CellModels that are active
    %       t1Transition - T1 transition
    %       merge_vertices_initializing - merges vertices closer than a
    %          certain threshold (use during constructing only)
    %       make_vertices - gets rid of vertices beyond image, as well
    %           as return which cell touches a Vertex (only when
    %           constructing)
    %   
    %   --- Visualize ---
    %       draw - draw single tissue configuration in image
    %       movie - returns a (Xs x Ys x T) image stack of tissue
    %              configuration as a function of time
    %
    % xies@mit.edu March 2015
    properties
        
        % -- current config --
        centroids % Ncx2 array of centroid locations
        vert_coords % Nvx2 array of vertex locations
        vertices % hashmap of Vertex.m objects
        cells % hashmap of the CellModels contained by the tissue
        interfaces % hashmap of Interface
        connectivity % adjacency matrix
        interVertDist % dist matrix
        parameters % simulation parameters
        contractile_params % contraction model + params
        
        Xs % tissue pixel size
        Ys
        merge_threshold_in_px
        t % timestamp
        
    end
    methods
        
        function tis = Tissue(regions,vert_coords,centroids,conn_opt,t)
            % Constructor for Tissue object
            %
            % USAGE
            %  (to create from scratch):
            %       tis = Tissue(regions, vert_coords, centroids,'purse string',t)
            %       tis = Tissue(regions, vert_coords, centroids,'purse string')
            %                       assumes t = 0
            %  (to copy object)
            %       tis = Tissue(old_tissue);
            %
            % March 2015, xies@mit.edu
            if nargin > 0
                if nargin < 5
                    t = 0;
                end
                
                if nargin == 1
                    % If there is a single input value, we must only want
                    % to copy an old tissue object (useful because
                    % container.Map is a reference object, and copies will
                    % modify the original)
                    
                    % Copy the value objects
                    tis_old = regions;
                    tis.centroids = tis_old.centroids;
                    tis.Xs = tis_old.Xs; tis.Ys = tis_old.Ys;
                    tis.merge_threshold_in_px = tis_old.merge_threshold_in_px;
                    tis.t = tis_old.t;
                    tis.vert_coords = tis_old.vert_coords;
                    tis.connectivity = tis_old.connectivity;
                    tis.interVertDist = tis_old.interVertDist;
                    tis.parameters = tis_old.parameters;
                    tis.contractile_params = tis_old.contractile_params;
                    
                    % Copy the reference object via concatenation w/ empty
                    % Map, creating a brand new Map!
                    % DO NOT REMOVE
                    tis.cells = [tis_old.cells; containers.Map()];
                    tis.vertices = [tis_old.vertices; containers.Map()];
                    tis.interfaces = [tis_old.interfaces; containers.Map()];
                    
                else % Else construct from scratch
                    
                    tis.t = t;
                    tis.centroids = centroids;
                    tis.Xs = size(regions,1); tis.Ys = size(regions,2);
                    
                    tis.vertices = containers.Map('KeyType','int32','ValueType','any');

                    % Merge vertices that are too close to each other
                    % @todo: This requires speedup
                    tis.merge_threshold_in_px = 3;
                    vert_coords = tis.merge_vertices_initializing( ...
                        vert_coords, tis.merge_threshold_in_px);
                    tis.vert_coords = vert_coords;
                    
                    % Instantiate valid vertices
                    tis = tis.make_vertices(regions);
                    
                    % Initialize a Map object to hold CellModel
                    num_cells = max(unique(regions));
                    tis.cells = containers.Map('KeyType','int32','ValueType','any');
                    
                    % Instatiate valid cells
                    vIDlist = tis.vertices.keys; vIDlist = [vIDlist{:}];
                    for i = 1:num_cells
                        tis.cells(i) = ...
                            CellModel(int32(i), tis,...
                            vIDlist( ...
                            cellfun(@(x) any(x == i), {tis.getVertices.cellIDs})), ... % Vertex keys
                            centroids(i,:) );
                    end
                    
                    % Make interfaces
                    tis = tis.connect_interfaces(conn_opt);
                    tis.interVertDist = squareform( pdist(tis.vert_coords) );
                    
                end
                
                if ~isValid(tis)
                    error('Tissue not self-consistent, exiting');
                end
                
            end
        end % Constructor
        
        function flag = isValid(tis)
            % Run self-consistency tests on current tissue state
            % Checks the following:
            %   1) vertices and vert_coords match and are ordered correctly
            %   2) vertices and the unique set of cell-owned vertices also
            %      match (not ordered)
            %   3) dimensions of connectivity and interVertDist are
            %      consistent with vert_coords
            
            flag = 1;
            % vertices and vert_coords match
            vx = [tis.getVertices.x]; vy = [tis.getVertices.y];
            flag = flag & all(all(cat(2,vx',vy') == tis.vert_coords));
            
            % Vertices and the set of cell vertices match
            c = tis.getCells;
            vt = tis.getVertices( unique([c.vIDs]) );
            vx = unique([vt.x]); vy = unique([vt.y]);
            flag = flag && all(vx == unique( [tis.getVertices.x] ) );
            flag = flag && all(vy == unique( [tis.getVertices.y] ) );
            
            if ~isempty( tis.connectivity)
                flag = flag && tis.vertices.length == mean(size(tis.connectivity));
                flag = flag && tis.vertices.length == mean(size(tis.interVertDist));
            end
            
        end
        
        % ------  Calculate energy, force, velocity ------
        
        function E = get_energy(tis)
            % GET_ENERGY Returns the current energy of the system
            % 
            % USAGE: E = get_energy(tis)
            %
            % Right now implements area elasticity, parameter elasticity,
            % and active contractility.
            
            
%             p = tis.parameters; % get parameters
% %             lambda = p.lengthScale;
            e = tis.interfaces.values; e = [e{:}];
            bondTerm = sum([e.energy]);
%             
            c = tis.cells.values(); c = [c{:}];
            cellTerm = sum([c.energy]);
            
            E = bondTerm + cellTerm;
%             
%             dist = tis.interVertDist .* tis.connectivity;
%             lineTensionTerm = nansum(nansum( triu(dist / lambda) ) * p.lineTension);
%             
%             current_areas = [tis.getCells.area];
%             areaElasticTerm = nansum( p.areaElasticity.* ...
%                 ([tis.getCells.area] - p.targetAreas).^2 / 2 /p.lambda^2 );
%             perimElasticTerm = nansum( p.perimElasticity.* ...
%                 ([tis.getCells.perimeter] - p.targetPerimeters) ./ p.lambda );
%             
%             C = [tis.getCells.contractility];
%             activeContractionTerm = nansum( C .* current_areas ./ p.lambda^2 );
%             
%             E_lat = areaElasticTerm + perimElasticTerm + ...
%                 lineTensionTerm;
%             E_cont = activeContractionTerm;
%             
%             E = E_lat + E_cont;

        end % get_energy
        
        function V = get_force(tis)
            % GET_FORCE
            %
            % Returns the forces as given by -grad(E) of current tissue
            % configuration.
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
            % The graduent terms are implemented as Nagai (2000).
            
            % Grab data
            vcoords = tis.vert_coords;
            D = tis.interVertDist;
            
            % Grap tissue-constant parameters
            length_scale = tis.parameters.lengthScale;
            
            % Initialize
            num_verts = size(vcoords,1);
            V = zeros( size(vcoords) );
            
            vIDList = tis.vertices.keys;
            vIDList = [vIDList{:}];
            % For now, loop through; vectorize later?
            for i = 1:num_verts
                
                % Find connected vertices
                vi = tis.vertices(vIDList(i));
%                 J = find(conn(i,:) == 1); % Idx of connected vertices
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
                    
%                     tis.draw('showVectors',{-line_tension_term,i},'showActive');
%                     keyboard
                    
                    % Go through all CELLS associated with current vertex,
                    % and calculate cell elasticity.
                    % NOTE that we can't just go through "edges" themselves
                    for this_cell = neighbors
                        
                        % Need to sort vertices counter-clockwise
                        sortedVt = tis.getVertices( this_cell.vIDs );
                        sortedVt = sortedVt.sort( this_cell.centroid );
                        
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
                        
                        % DEBUGGING
%                         tis = tis.activateCell( this_cell.cellID );
%                         tis.draw('showVectors',{v,i},'showActive');
%                         tis = tis.deactivateCell( this_cell.cellID );
%                         hold on, sortedVt(I).draw;
                        
                    end
                    
                    V(i,:) = line_tension_term + area_elastic_term + ...
                        perim_elastic_term + active_contraction_term;
                    
%                     tis.draw('showVectors',{V(i,:),i},'showActive');
%                     drawnow

                end
                
                if any(any(isnan( V ))), keyboard; end
                
            end
            
            % Check that fixed vertices have not moved
            if any( V(tis.parameters.fixed_verts,:) > 0 )
                keyboard;
            end
            
        end % get_velocities
        
        % ------ Measurements (area, contractility, etc) ---------
        
        function A = getArea(tisArray)
            % Returns an array (or vector) of the cell areas contained in
            % input tissue array.
            T = numel(tisArray);
            % @todo: assumes tisArray does not change number of cells over
            % time!!
            % Idea for imporving: use cellID to make a NaN-padded array?
            num_cells = tisArray(1).cells.length;
            A = nan( T,num_cells );
            for i = 1:T
                c = tisArray(i).getCells;
                c = c.sortbyID;
                A(i,:) = [c.area];
            end
            A = A*tisArray(1).parameters.um_per_px^2;
        end % getArea
        
        function C = getContractility(tisArray)
            % Returns an array (or vector) of the contractility contained in
            % input tissue array.
            T = numel(tisArray);
            % @todo: assumes tisArray does not change number of cells over
            % time!!
            % Idea for imporving: use cellID to make a NaN-padded array?
            num_cells = tisArray(1).cells.length;
            C = nan( T,num_cells );
            for i = 1:T
                c = tisArray(i).getCells;
                c = c.sortbyID;
                C(i,:) = [c.contractility];
            end
        end % getContractility
        
        function CX = getCentroidX(tisArray)
            % Returns an array (or vector) of the centroid-x contained in
            % input tissue array.
            T = numel(tisArray);
            % @todo: assumes tisArray does not change number of cells over
            % time!!
            % Idea for imporving: use cellID to make a NaN-padded array?
            num_cells = tisArray(1).cells.length;
            CX = nan( T,num_cells );
            for i = 1:T
                c = tisArray(i).getCells;
                c = c.sortbyID;
                cx = cat(1,c.centroid);
                CX(i,:) = cx(:,2);
            end
        end % getCentroidX
        
        function CY = getCentroidY(tisArray)
            % Returns an array (or vector) of the centroid-y contained in
            % input tissue array.
            T = numel(tisArray);
            % @todo: assumes tisArray does not change number of cells over
            % time!!
            % Idea for imporving: use cellID to make a NaN-padded array?
            num_cells = tisArray(1).cells.length;
            CY = nan( T,num_cells );
            for i = 1:T
                c = tisArray(i).getCells;
                c = c.sortbyID;
                cy = cat(1,c.centroid);
                CY(i,:) = cy(:,1);
            end
        end % getCentroidY
        
        function T = getTime(tisArray)
            % Returns a vector of the simulation time contained in
            % input tissue array.
            T = [tisArray.t];
        end % getTime
        
        % ------ Simulation bookkeeping methods ---------
        
        function tis = evolve(tis_old, new_vcoords, varargin)
            % EVOLVE
            % Updates and returns a new copy of the old tissue
            % configuration by moving all the vertex positions
            %   NOTA BENE: the old .cells container is COPIED and not
            %              direclty modified, since it's a reference object
            %
            % USAGE: new_tissue = old_tissue( new_vcoords );
            %        new_tissue = old_tissue( new_vcoords ,'no_update');
            
            % RESHAPE new_coords if fed in from ODE solver
            if isvector(new_vcoords)
                N = numel(new_vcoords);
                new_vcoords = reshape(new_vcoords,N/2,2);
            end
            
            if size(new_vcoords,1) ~= tis_old.vertices.length
                error('Size of new vertex list must match old vertices')
            end
            tis = Tissue(tis_old);
            tis.vert_coords = new_vcoords;
            vIDList = tis.vertices.keys; vIDList = [vIDList{:}];
            
            for i = 1:numel(vIDList)
                % Move vertices in CellModels
                v = tis.vertices(vIDList(i));
                cContainV = v.cellIDs;
                for c = cContainV'
                    tis.cells( c ) = ...
                        tis.cells(c).updateCell( tis );
                end
                
                % Update edges associated with this vertex
                eContainsV = v.bondIDs;
                for e = ensure_row(eContainsV)
                    tis.interfaces( e) = ...
                        tis.interfaces(e).updateInterface(tis);
                end
                
                % Move Vertex
                tis.vertices(vIDList(i)) = ...
                    tis.vertices(vIDList(i)).move(new_vcoords(i,:));
            end
            
            % Update distance maps
            D = squareform(pdist(tis.vert_coords));
            tis.interVertDist = D;
            
%             % Perform T1 transitions
            D(~tis.connectivity) = NaN;
            D(tis.parameters.fixed_verts,:) = NaN;
            D = triu(D);
            D( D==0 ) = NaN;
            [I,J] = find( D < tis.parameters.t1Threshold);
            for i = 1:numel(I)
                tis = tis.t1Transition( tis.getVertices(vIDList([I(i) J(i)])) );
            end
            
            % Consistency check
            if ~tis.isValid, keyboard; end
            
            % Advance time stamp by one
            if nargin < 3,
                tis.t = tis.t + tis.parameters.dt_per_frame;
            end
            
        end % evolve
        
        function tis = updateVertCoords(tis)
            % If you move a vertex without moving vert_coords explicitly,
            % use this to make them consistent.
            v = tis.getVertices;
            vx = [v.x]; vy = [v.y];
            tis.vert_coords = cat(2,vx',vy');
        end
        
        function tis = setParameters(tis,varargin)
            % Sets the simluation/evolution parameters
            %
            % USAGE: tissue = 
            %           tis.setParameters('parameter',value);
            %   PARAMETERS:
            %       'lineTension' - scalar or Nv x Nv matrix, no defaults
            %       'lineAnisotropy' - anisotropic tension vertical v.
            %                      horizontal junctions, default = 1
            %       'areaElasticity' - scalar or Nc x 1 vector, no
            %                       defaults
            %       'perimeterElasticity' - scalar or Nc x 1 vector, no
            %                       defaults
            %       'connectivity' - string: see Matrix
            %       'targetArea' - scalar or Nc x 1 vector, default is the
            %                       mean value of current setting
            %       'targetPerimeter' - default = 0
            %       'dimensonless' - 1 for reduced parameter values,
            %                       default = 0
            %       'lengthScale' - for reducing parameters
            %       'forceScale' - usually equal to linetension unless
            %              there is special cases set or anisotropy
            %       'viscosity' - scalar, no defaults
            %       -- other parameters
            %       'jitterSize' - for jittering vertex locations (px)
            %                       default = 0
            %       't1Threshold' -- default = 1 (px!)
            %       'timeStep' - for Euler integration, default = 0.001
            %       'frame_per_sec' - default = 1
            %       'um_per_px' - default = 1
            % 
            % xies@mit April 2015
            
            p = myParameterParser;
            Nb = tis.interfaces.length;
            Nc = tis.cells.length;
            
            % --- Dimensonality ---
            % NO IDEA why matlab throws errors all over the place if you
            % have a validator here...
            addOptional(p,'dimensionless',false,@(x) true); %default = 0
            
            % --- length scale ---
            addOptional(p,'lengthScale',1,@isscalar); %default = 1
            
            % --- force scale ---
            addOptional(p,'forceScale',1,@isscalar); %default = 1
            
            % --- Area elasticity ---
            addRequired(p,'areaElasticity', ...
                @(x) isscalar(x) || (isvector(x) && numel(x) == Nc) );
            
            % --- Line tension ---
            addRequired(p,'lineTension',...
                @(x) isscalar(x) || (isvector(x) && numel(x) == Nb) );
            
            % --- Line anisotropy ---
            addOptional(p,'lineAnisotropy',1,@isscalar);
            
            % --- Perimeter elasticity ---
            addRequired(p,'perimElasticity', ...
                @(x) isscalar(x) || (isvector(x) && numel(x) == Nc) );
            
            % --- Viscosity ---
            addRequired(p,'viscosity',@isscalar);
            
            % --- Connectivity ---
            validConnectivityOpts = {'purse string'};
            addRequired(p,'connectivity',...
                @(x) any(validatestring(x,validConnectivityOpts)));
            
            % --- target area ---
            c = tis.getCells; % default = mean of current areas
            addOptional(p,'targetArea',mean([c.area]), ...
                @(x) isscalar(x) || (isvector(x) && numel(x) == Nc) );
            
            % --- target perimeter ---
            addOptional(p,'targetPerimeter',0, ... % default = 0
                @(x) isscalar(x) || (isvector(x) && numel(x) == Nc) );
            
            % --- bookkeeping ---
            addOptional(p,'jitterSize',0,@isscalar); %default = 0
            addOptional(p,'stepSize',0.001,@isscalar); %default = 0.001
            addOptional(p,'t1Threshold',1,@isscalar) % default = 1
            
            % --- Units ---
            addOptional(p,'um_per_px',1,@isscalar); %default = 1
            addOptional(p,'dt_per_frame',1,@isscalar); %default = 1
            
            parse(p,varargin{:});
            
            tis.parameters = p.Results;
            tis.parameters.fixed_verts = cellfun( ...
                @numel,{tis.getVertices.cellIDs}) < 3;
            
            if tis.parameters.dimensionless
                % Find dimensional parameters
                sigma = tis.parameters.forceScale; % force dimension
                lambda = tis.parameters.lengthScale; % length simension
                K_a = tis.parameters.areaElasticity * sigma / lambda^3;
                K_p = tis.parameters.perimElasticity * sigma / lambda;
                A0 = lambda^2 * 3 * sqrt(3) / 2;
            else
                sigma = tis.parameters.lineTension;
                K_a = tis.parameters.areaElasticity;
                K_p = tis.parameters.perimElasticity;
                A0 = tis.parameters.targetAreas;
            end
            
            % Set Intarface-related parameters
            tis = tis.setLineTension;
            if tis.parameters.lineAnisotropy ~= 1
                tis = tis.setLineAnisotropy;
            end
            
            % Set CellModel parameters
            tis = tis.setAreaElasticity;
            tis = tis.setPerimElasticity;
            tis = tis.setTargetArea;
            
            display('---Parameter check---')
            display(['Normalized line tension = '...
                num2str( sigma / K_a / A0^3/2)])
            display('(Should be ~0.12)')
            display('-')
            display(['Normalized perimeter tension = '...
                num2str( K_p / K_a / A0)])
            display('(Should be ~0.04)')
            display('---------------------')
            
        end % setParameters
        
        function tis = activateCell( tis, cellIDs, varargin)
            % activateCell
            % Set specified cells (IDs) to "active = 1". Used to keep track
            % of who's "ventral" in the model.
            % Usage: tis = tis.activateCell(cellIDs)
            %        tis = tis.activateCell(cellIDs, alt_tension)
            % Optionally, set non-active cells' lineTension to a specified
            % value. (See Spahn 2014).
            
            % If no arguments are given, then activate all cells
            if nargin < 2,
                cellIDs = tis.cells.keys;
                cellIDs = [cellIDs{:}];
            end
            % Activate specified cells
            for i = 1:numel(cellIDs)
                tis.cells( cellIDs(i) ) = ...
                    tis.cells( cellIDs(i) ).activateCell;
            end
            
            if nargin > 2 % Specified alt tension
                nonActiveTension = varargin{1};
                inactiveCells = tis.getInactiveCells;
                num_cells = numel(inactiveCells);
                for i = 1:num_cells
                    
                    edges = tis.getInterfaces( inactiveCells(i).bondIDs );
                    tis = tis.setLineTension( [edges.ID], nonActiveTension );
                    
                end
            end
            
        end % activateCell
        
        function tis = deactivateCell(tis, cellIDs)
            % Deactivates cell(s)
            % USAGE: tis = tis.deactivateCell; (default = all cells)
            %        tis = tis.deactivateCell(IDs);
            if nargin < 2 % If no cells specified, deactivate all cells
                cellIDs = tis.cells.keys;
                cellIDs = [cellIDs{:}];
            end
            for i = 1:numel(cellIDs)
                tis.cells( cellIDs(i) ) = ...
                    tis.cells( cellIDs(i) ).deactivateCell;
            end
        end % deactivateCell
        
        function tis = deactivateBorder(tis,n)
            % deactivateBorder:
            % Deactivates the border (up to n order) cells
            % USAGE: tis = tis.deactivateBorder(n)
            %
            % @todo: implement orders n > 1
            
            % By default only deactivate first layer
            if nargin < 2, n = 1; end
            
            if n>1, error('Haven''t implemented higher order :('); end
            
            cellIDList = tis.cells.keys; cellIDList = [cellIDList{:}];
            for ID = cellIDList
                if tis.numCellNeighbors( tis.cells(ID) ) < 6
                    tis = tis.deactivateCell( ID );
                end
            end
            
        end % deactivateBorder
        
        function tis = setContractility(tis,C,cIDs)
            % Directly sets the active contractility coefficient
            % Requires all cells to have its own specified contractility
            % 
            % USAGE: tis = tis.setContractility( C , cIDs );
            %        tis = tis.setContractility( C , cIDs );
            % INPUT: tis - tissue
            %        C - vector of contractility values
            %        cIDs - by default sets all cells
            
            if nargin < 3
                cIDs = tis.cells.keys; cIDs = [cIDs{:}];
            end
            
            if numel(cIDs) ~= numel( C )
                error('Number of contractility coeff and number of cells don''t match');
            end
            
            for i = 1:numel(cIDs)
                if tis.cells(cIDs(i)).isActive
                    % set cell contractility iff it's an active cell
                    tis.cells( cIDs(i) ) = ...
                        tis.cells( cIDs(i) ).setContractility(C(i));
                end
            end
        end % setContractility
        
        function tis = setContractilityModel(tis,modelfun,params,cIDs)
            % setContractilityModel - sets the value of contractility
            % according to the given model.
            % 
            % USAGE: tis = tis.setContractilityModel( modelfun, params,cIDs)
            %        tis = tis.setContractilityModel( modelfun, params)
            %                   - sets contractility in all ACTIVE cells
            %
            % INPUT: modelfun - a function handle, should return the
            %                   contractility if given the centroid
            %                   locations (Nc x 2) and the parameters also
            %                   passed to this function.
            %                   
            %        params
            %
            % Records the modelfun + params in tis.contraction_param
            
            if nargin < 4
                cellsOI = tis.getActiveCells;
            else
                cellsOI = tis.getCells( cIDs );
            end
            
            % Evaluate modelfun and set the values
            ct = cat(1,cellsOI.centroid);
            C = modelfun(ct,params);
            tis = tis.setContractility( C,[cellsOI.cellID] );
            tis = tis.deactivateBorder;
            
            % Record modelfun + param
            tis.contractile_params.model = modelfun;
            tis.contractile_params.params = params;
            
        end
        
        function tis = jitterVertices(tis, STD)
            % Add a set amount of Gaussian jitter to vertex position.
            %
            % USAGE: tis = tis.jitterVertices( STD )
            
            verts = tis.vert_coords;
            I = ~tis.parameters.fixed_verts;
            jitter = STD*randn([ numel(I(I)), 2]);
            verts( I,: ) = verts( I,: ) + jitter;
            
            tis = tis.evolve( verts , 'no_update' );
            
        end
        
        % ------ Verted-vertex connectivity ------
        
        function tis = connect_interfaces(tis, opt, cells)
            % Sets the vertex-vertex connectivity matrix according to the
            % specified model configurations, and then makes the
            % appropriate Interface objets
            % 
            % USAGE: 
            % tis = tis.connect_interfaces( opt )
            % tis = tis.connect_interfaces( opt )
            % 
            % INPUT: tis - the tissue to be connected
            %        opt - 'purse string' / 'apical' / 'both'
            % 
            % @todo: implement 'apical' and 'both'
            
            vIDsList = tis.vertices.keys; vIDsList = [vIDsList{:}];
            num_vertices = numel(vIDsList);
            % Start with brand-new hashmap
            tis.interfaces = containers.Map('KeyType','int32','ValueType','any');
            
            switch opt
                % Connect the 'junctions' of cells only
                case 'purse string'
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
                                % Check for double counting:
                                J = find( vIDsList == connectedVertIDs(j) );
                                if conn(i,J) || conn(J,i), continue; end
                                
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
                    
                otherwise
                    error('Unrecognized vertex connection option.')
            end
            
            tis.connectivity = conn;
            
        end % connectVertices
        
        % ------ Cell-Vertex connectivity -----
        
        function cellsThatTouch = cellsContainingVertex(tis,vert)
            % Return a list of CellModel that contains the current Vertex
            %
            % Usage: touchingCells = tis.cellsContainingVertex( vert )
            
            if numel(vert) > 1, error('Can''t do more than 1 vertex.'); end
            N = numel(vert.cellIDs);
            cellsThatTouch(1:N) = CellModel();
            for i = 1:N
                cellsThatTouch(i) = tis.cells( vert.cellIDs(i) );
            end
            
            if isempty(cellsThatTouch), keyboard; end
        end % cellsContainingVertex
        
        % ----- Cell-cell connectivity ----
        
        function neighbors = neighborsOfCell( tis, cellOI, order_n )
            % Return a annulus of cells around cellOI, at order_n-th layer
            % USAGE: neighbors = neighborsOfCell( tis, cellOI, order_n )
            
            % If input is only a cellID instead of object itself, then fetch cell
            if isnumeric(cellOI), cellOI = tis.cells(cellOI); end
            
            tmpCellSet = java.util.HashSet;
            % Add the n-th neighbors
            nth_neighb = tis.allNthOrderNeighbors(cellOI, order_n);
            for c = nth_neighb
                tmpCellSet.add( c.cellID );
            end
            % Delete the n-1-th neighbors
            nminus1_neighbor = tis.allNthOrderNeighbors(cellOI, order_n-1);
            for c = nminus1_neighbor
                tmpCellSet.remove( c.cellID );
            end
            
            % Grab neigbor cells from the Map in Tissue
            tmpCellSet = tmpCellSet.toArray;
            neighbors(1:numel(tmpCellSet)) = CellModel;
            for i = 1:numel(tmpCellSet)
                neighbors(i) = tis.cells(tmpCellSet(i));
            end
            
        end %neighborsOfCell
        
        function neighbors = allNthOrderNeighbors( tis, cellOI, order_n )
            % allNthOrderNeighbors - Find the "halo" of cells around cellOI
            % until we've filled it to ORDER_N layers.
            %
            % USAGE: neighbors = allNthOrderNeighbors( tis, cellOI, orderN )
            %
            % OUTPUT: neighbors - cellModel array
            %
            % Not very fast -- @todo: need HashSet to work with matlab custom
            % classes, can't figure out yet.
            
            % If input is only a cellID instead of object itself, then fetch cell
            if isnumeric(cellOI), cellOI = tis.cells(cellOI); end
            
            tmpCellSet = java.util.HashSet; % HashSet constant add time
            % Add the central cell; can only use primitive data types with
            % imported HashSet :(
            cellIDList = tis.cells.keys();
            tmpCellSet.add(cellOI.cellID);
            for i = 1:order_n
                
                % Grab all the current i-th order filled-in corona
                tmpCells = tmpCellSet.toArray;
                for d = 1:numel(tmpCells)
                    % Iterate through all cells to find all connected cells
                    thisCell = tis.cells( tmpCells(d) );
                    %                     if i == 2, keyboard; end
                    for c = 1:tis.cells.length
                        if tis.connected( tis.cells(cellIDList{c}), thisCell )
                            tmpCellSet.add( tis.cells(cellIDList{c}).cellID );
                        end
                    end
                end
                
            end
            
            % Grab neigbor cells from the Map in Tissue
            tmpCellSet = tmpCellSet.toArray;
            neighbors(1:numel(tmpCellSet)) = CellModel;
            for i = 1:numel(tmpCellSet)
                neighbors(i) = tis.cells(tmpCellSet(i));
            end
            
        end % allNthOrderNeighbors
        
        function numNeigh = numCellNeighbors(tis, cellA)
            % Return the number of neighboring cells to input cell
            % USAGE: numNeighbor = tis.numCellNeighbors(input_cell)
            numNeigh = numel( tis.neighborsOfCell( cellA,1 ) );
        end % numCellNeighbors
        
        function flag = connected(~,cell_a, cell_b)
            % Basic function where if two cells share one vertex,
            % they're connected
            % 
            % Usage: flag = tis.connected( cells_a, cell_b)
            
            flag = 0;
            va = cell_a.vIDs;
            for i = 1:numel(va)
                % compare vID / vertices key values
                flag = any(va(i) == cell_b.vIDs);
                if flag, return; end
            end
        end % connected
        
        % ----- Cell handling and measurements ------
        
        function cells = getCells( tis, varargin )
            % Returns cells from tissue as an array.
            % 
            % USAGE:
            %    cells = tis.getCells(); % returns all cells
            %    cells = tis.getCells(cellID); returns subset
            
            % Get all cells
            if nargin == 1
                cells = tis.cells.values;
                cells = [cells{:}];
            else
                cellID = varargin{1};
                num_cells = numel(cellID);
                cells(1:num_cells) = CellModel();
                for i = 1:num_cells
                    cells(i) = tis.cells( cellID(i) );
                end
            end
        end % getCells
        
        function cells = getActiveCells(tis)
            % Return all active cells in the tissue
            %
            % USAGE: actives = tissue.getActiveCells;
            cells = tis.getCells;
            cells = cells( [cells.isActive] > 0 );
        end % getActiveCells
        
        function cells = getInactiveCells(tis)
            % Return all non-active cells in the tissue
            %
            % USAGE: actives = tissue.getInactiveCells;
            cells = tis.getCells;
            cells = cells( [cells.isActive] == 0 );
        end % getActiveCells
        
        function cIDs = getCellsWithinRegion(tis, box)
            % Returns all cells whose centroid lies within a region
            % specified by BOX = [X0 Y0 Xf Yf] (inclusive)
            %
            % USAGE: cellIDs = tis.getCellsWithinRegion( box )
            
            cells = tis.getCells;
            ct = cat(1,cells.centroid);
            cellIDList = tis.cells.keys();
            
            I = ct(:,1) >= box(1) & ct(:,1) <= box(3);
            I = I & ct(:,2) >= box(2) & ct(:,2) <= box(4);
            
            cIDs = [cellIDList{I}];
            
        end % getCellsWithinRegion
        
        function tis = setTargetArea(tis,cIDs,alt_area)
            % Passes parameters from Tissue into CellModels
            % 
            % Usage: tis = setTargetArea(tis,cIDs,alt_area)
            %
            % INPUT:
            %   cIDs - by default does every cell
            %   alt_area - should match size of cIDs
            
            if nargin < 2
                cIDs = tis.cells.keys;
                cIDs = [cIDs{:}];
                alt_area = tis.parameters.targetArea;
            end
            
            num_cells = numel(cIDs);
            % Update targetAreas
            for i = 1:num_cells
                c = tis.cells( cIDs(i) );
                if isscalar(alt_area)
                    c.targetArea = alt_area;
                else
                    c.targetArea = alt_area(i);
                end
                tis.cells( cIDs(i) ) = c;
            end
        end
        
        function tis = setAreaElasticity(tis,cIDs,alt_ka)
            % Passes parameters from Tissue into CellModels
            % 
            % Usage: tis = setAreaElasticity(tis,cIDs,alt_ka)
            %
            % INPUT:
            %   cIDs - by default does every cell
            %   alt_area - should match size of cIDs
            
            if nargin < 2
                cIDs = tis.cells.keys;
                cIDs = [cIDs{:}];
                alt_ka = tis.parameters.areaElasticity;
            end
            
            num_cells = numel(cIDs);
            % Update targetAreas
            for i = 1:num_cells
                c = tis.cells( cIDs(i) );
                if isscalar(alt_ka)
                    c.areaElasticity = alt_ka;
                else
                    c.areaElasticity = alt_ka(i);
                end
                tis.cells( cIDs(i) ) = c;
            end
        end
        
        function tis = setPerimElasticity(tis,cIDs,alt_kp)
            % Passes parameters from Tissue into CellModels
            % 
            % Usage: tis = setPerimElasticity(tis,cIDs,alt_kp)
            %
            % INPUT:
            %   cIDs - by default does every cell
            %   alt_kp - should match size of cIDs
            
            if nargin < 2
                cIDs = tis.cells.keys;
                cIDs = [cIDs{:}];
                alt_kp = tis.parameters.perimElasticity;
            end
            
            num_cells = numel(cIDs);
            % Update targetAreas
            for i = 1:num_cells
                c = tis.cells( cIDs(i) );
                if isscalar(alt_kp)
                    c.perimElasticity = alt_kp;
                else
                    c.perimElasticity = alt_kp(i);
                end
                tis.cells( cIDs(i) ) = c;
            end
        end
        
        % ----- Interface handling ----
        
        function e = getInterfaces(tis, bIDs)
            % Returns interfaces from the tissue with specified IDs
            % Else returns all.
            % 
            % USAGE: e = tis.getInterfaces();
            % USAGE: e = tis.getInterfaces(bIDs);
            
            if nargin < 2,
                bIDs = tis.interfaces.keys;
                bIDs = [bIDs{:}];
            end
            
            num_edges = numel(bIDs);
            e(1:num_edges) = Interface;
            for i = 1:num_edges
                e(i) = tis.interfaces( bIDs(i) );
            end
        end
        
        function tis = setLineTension(tis,bIDList,L)
            % Set the line tension of each bond in the tissue.
            % Can specify each bond or have them be all the same
            %
            % USAGE: tis = tis.setLineTension(bIDs,L)
            %        tis = tis.setLineTension - uses tis.parameters
            %
            % INPUT:
            %   tis.parameters.lineTension:
            %       scalar or Nv by Nv matrix with the same sparsity as
            %       tis.connectivity
            
            if nargin < 3
                bIDList = tis.interfaces.keys; bIDList = [bIDList{:}];
                L = tis.parameters.lineTension;
            end
            
            for i = 1:numel(bIDList)
                e = tis.interfaces( bIDList(i) );
                if isscalar(L)
                    e.tension = L;
                else
                    e.tension = L(i);
                end
                tis.interfaces( bIDList(i) ) = e;
            end
            if nargin > 1
                allIDs = tis.interfaces.keys; allIDs = [allIDs{:}];
                tis.parameters.lineTension( ismember(bIDList,allIDs) ) = L;
            end
            
        end % setLineTension
        
        function tis = setLineAnisotropy(tis,bIDList,a)
            % Set the line tension anisotropically according to
            % orientation.
            %
            % Anisotropy value is the same across specified cells
            %
            % USAGE: tis = tis.setLineAnisotropy(bIDs,a)
            %        tis = tis.setLineAnisotropy - use tis.parameters as
            %           default value
            
            sigma = tis.parameters.forceScale;
            if nargin < 3
                a = tis.parameters.lineAnisotropy;
                bIDList = tis.interfaces.keys; bIDList = [bIDList{:}];
            end
            
            if a == 1, return; end % Do nothing if a is unit
            
            new_tensions = zeros(1,numel(bIDList));
            for i = 1:numel(bIDList)
                e = tis.interfaces( bIDList(i) );
                theta = e.angle;
                % Scale horizontal junctions by factor
                if (theta > 0 && theta < pi/4) ...
                        || (theta > 3*pi/4 && theta < pi)
                    e.tension = sigma*a;
                end
                % set new junctions
                tis.interfaces( bIDList(i) ) = e;
                new_tensions(i) = e.tension;
            end
            
            % Update parameter list
            allIDs = tis.interfaces.keys; allIDs = [allIDs{:}];
            if isscalar( tis.parameters.lineTension )
                tis.parameters.lineTension = ...
                    tis.parameters.lineTension*ones(1,numel(allIDs));
            end
            tis.parameters.lineTension( ismember(bIDList,allIDs) ) = new_tensions;
            
        end % setLineAnisotropy
        
        % ----- Vertex handling ------
        
        function tis = t1Transition( tis, vt )
            % Make a T1 transition (a la Lecuit definition)
%             figure, tis.draw('showVertices');
            if numel(vt) > 2
                keyboard
            end
            
            % Swap cell ownership
            % figure out which cells contain each vertex
            cells1 = vt(1).cellIDs;
            cells2 = vt(2).cellIDs;
            % Figure out which bonds contain both vertices (should be only
            % one!)
            edge = tis.interfaces( intersect(vt(1).bondIDs, vt(2).bondIDs) );
            
            % figure out which cells contain both or just one
            cells_both = tis.getCells( intersect(cells1,cells2) );
            cells_single = ...
                tis.getCells( setxor(cells1,cells2) );
            
            if numel(cells_both) > 2 || numel(cells_single) > 2
                keyboard % ill-defined!
            end
            
            % Rotate vertex positions and add both cells_single IDs
            midpt(1) = mean([vt.x]); midpt(2) = mean([vt.y]);
            for i = 1:numel(vt)
                vt(i) = vt(i).rotate(midpt,-pi/2);
                tis.vertices( vt(i).ID ) = vt(i);
                vt(i).cellIDs = union( vt(i).cellIDs, [cells_single.cellID]);
            end
            tis = tis.updateVertCoords;
            
            % Swap interface cell ownership and update tissue
            edge.cIDs = [cells_single.cellID];
            tis.interfaces( edge.ID) = edge.updateInterface(tis);
            
            % Cells with only 1 before now has 2 and the bond
            for i = 1:numel(cells_single)
                % Add vert to cell list
                cells_single(i).vIDs = union( cells_single(i).vIDs, [vt.ID]);
                % Add bond to cell list
                cells_single(i).bondIDs = union( ...
                    cells_single(i).bondIDs, edge.ID );
                tis.cells( cells_single(i).cellID ) = ...
                    cells_single(i).updateCell(tis);
            end
            
            % Figure out which vertex is closer to each cells_both
            for i = 1:numel(cells_single)
                % Find the vertex that's farther away
                D = cells_both(i).get_distance_to( [[vt.x];[vt.y]] );
                [~,I] = max(D);
                % remove the vID from cellModel
                cells_both(i).vIDs = setdiff(cells_both(i).vIDs, vt(I).ID );
                % remove bondID from cellModel
                cells_both(i).bondIDs = setdiff(cells_both(i).bondIDs, edge.ID );
                % remove cellID from vertex Model
                vt(I).cellIDs = setdiff( vt(I).cellIDs, cells_both(i).cellID );
                % Put vertex and cells back in tissue
                tis.vertices( vt(I).ID ) = vt(I);
                tis.cells( cells_both(i).cellID ) = ...
                    cells_both(i).updateCell(tis);
            end
            
            % Update matrices
            tis.interVertDist = squareform(pdist( tis.vert_coords ));
            
            keyboard
        end
        
        function tis = make_vertices(tis, regions)
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
        
        function vcoords = merge_vertices_initializing(~,vcoords, ...
                merge_threshold_in_px)
            % Merge vertices that are closer than the specified threshold
            %
            % USAGE:
            % vcoords = tis.merge_vertices_initializing( vcoords,threshold)
            %
            % Use only from constructor
            
            num_vertices = size(vcoords);
            % Merge vertices which are super close to each other
            vertDist = squareform(pdist(vcoords));
            vertDist( logical(eye(num_vertices)) ) = NaN;
            
			% Loop merge function until no pairs of vertices are closer than
			% the threshold
            while any(any(vertDist <= merge_threshold_in_px))
                
                % Find a set of vertices to merge
                [I,J] = find(vertDist <= merge_threshold_in_px,1,'first');
                tobeMergedInd = [I,find(vertDist(I,:) <= merge_threshold_in_px)];
				mergedV = mean( vcoords(tobeMergedInd,:) );
                
                % Delete old ummerged vertices
                vcoords(tobeMergedInd,:) = [];
                % Add new merged vertex
                vcoords = cat(1,vcoords,mergedV );
                
                % update the distance maps
                num_vertices = size(vcoords,1);
                vertDist = squareform(pdist(vcoords));
                vertDist( logical(eye(num_vertices)) ) = NaN;
                
            end
        end % merge_vertices_initializing
        
        function verts = getVertices( tis , IDs )
            % Returns all or specified subset of vertices as array
            % USAGE: vts = tis.getVertices
            % 	     vts = tis.getVertices( IDs )
            
            % Return all Vertex objects as an array
            if nargin < 2
                verts = tis.vertices.values;
                verts = [verts{:}];
            else
                verts(1:numel(IDs)) = Vertex();
                for i = 1:numel(IDs)
                    verts(i) = tis.vertices( IDs(i) );
                end
            end
        end % getVertices
        
        %------ Visualization -----
        
        function I = draw(tis,varargin)
            % Draws a single tissue in binary image. Can return just the
            % outlines of cells, or also shade-in the active cells.
            %
            % USAGE: I = tis.draw;
            %        I = tis.draw('showActive');
            %           Highlights currentlly "active" cells
            %        I = tis.draw('showContractile')
            %           Highlights contractile cells and false-colors
            %           contractility amount.
            %        I = tis.draw('showVectors',V);
            %           Shows vector field V over all vertices
            %        I = tis.draw('showVectors',{v,ID});
            %           Shows only single vector v at vertex == ID
            %
            % Cannot handle more than one tissue. Use tisArray.movie() for
            % making a movie.
            
            if numel(tis) > 1, error('Can only handle single tissue; use tis.movie() to show movie.'); end
            
            I = zeros(tis.Xs,tis.Ys);
            cellIDList = tis.cells.keys();
            for i = 1:numel(cellIDList)
                I = I + tis.cells(cellIDList{i}).draw(tis);
                I = logical(I);
            end
            
            I = double(I) * 255;
            % Highlight active cells
            if any(strcmpi(varargin, 'showActive'))
                % Show active cells as filled-ins
                M = zeros(tis.Xs,tis.Ys);
                Acells = tis.getActiveCells;
                for i = 1:numel(Acells)
                    M = M + Acells(i).drawMask(tis);
                end
                M = M * 50;
                I = I + M;
            end
            
            % Highlight active cells
            if any(strcmpi(varargin, 'showContractile'))
                % Show active cells as filled-ins, with value proportional
                % to contractility
                M = zeros(tis.Xs,tis.Ys);
                Acells = tis.getActiveCells;
                C = [Acells.contractility];
                cmax = max(C);
%                 cmax = 0.05;
                for i = 1:numel(Acells)
                    M = M + Acells(i).drawMask(tis) * ...
                        Acells(i).contractility * 150 / cmax;
                end
                I = I + M;
            end
            
            % Highlight specific cells
            ind = find( strcmpi(varargin,'showCellID') );
            if ~isempty(ind)
                M = zeros(tis.Xs,tis.Ys);
                IDs = varargin{ind+1};
                Acells = tis.getCells(IDs);
                for i = 1:numel(Acells)
                    M = M + Acells(i).drawMask(tis);
                end
                I = I + M * 50;
            end
            
            hold off, imagesc(I), axis equal;
            
            % Highlight vertices
            if any(strcmpi(varargin, 'showVertices'))
                hold on
                vcoords = tis.vert_coords;
                scatter(vcoords(:,2),vcoords(:,1),'w');
            end
            
            % Display vector field
            ind = find( strcmpi(varargin,'showVectors') );
            if ~isempty(ind)
                if numel(varargin) < ind + 1, error('Need vector to draw'); end
                V = varargin{ind+1};
                % If input is not a cell object, then display all vectors
                if ~iscell(V)
                    hold on;
                    quiver(tis.vert_coords(:,2),tis.vert_coords(:,1), ...
                        V(:,2),V(:,1),0,'w-');
                else
                    % If it's a cell obj, then only the given vertex will
                    % have a vector over it
                    v = V{1};
                    ID = V{2};
                    if numel(ID) ~= size(v,1);
                        error('# of vectors should equal # of origins')
                    end
                    hold on
                    quiver(tis.vert_coords(ID,2),tis.vert_coords(ID,1), ...
                        v(:,2),v(:,1),0,'w-');
                end
            end
            
            hold off
            
        end
        
        function F = movie(tissues,varargin)
            % Make a movie of tissue evolving. Can return just the
            % outlines of cells, or also shade-in the active cells.
            %
            % USAGE: I = movie(tisues);
            %        I = movie(tisues,'showActive');
            
            num_frames = numel(tissues);
            if nargin == 1, opt = 'none';
            else
                opt = varargin{1};
            end
            
            F = zeros(tissues(1).Xs, tissues(1).Ys, num_frames);
            
            for f = 1:num_frames
                F(:,:,f) = tissues(f).draw(opt);
            end
            
        end
        
    end % Methods
    
    
end
