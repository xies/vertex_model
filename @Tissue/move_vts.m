function move_vts(tis, new_vcoords, new_time, varargin)
%Tissue/move_vts
% Updates and returns a new copy of the old tissue  configuration by
% moving all the vertex positions, updating Vertex, Interface, and
% CellModel properties, and then perform T1 transitions. Records current
% time and energy
%
% USAGE: tissue.move_vts( new_vcoords );
%        tissue.move_vts( new_vcoords ,'no_update');
% 

% RESHAPE new_coords if called from an ODE solver
if isvector(new_vcoords)
    N = numel(new_vcoords);
    new_vcoords = reshape(new_vcoords,N/2,2);
end

if size(new_vcoords,1) ~= tis.vertices.length
    error('Size of new vertex list must match old vertices')
end

tis.vert_coords = new_vcoords;
vIDList = tis.vertices.keys; vIDList = [vIDList{:}];

% Individually move vertices using Vertex.move
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

% Get new time stamp
tis.t = new_time;

% Advance T1 transition timestamp, delete any entries that
% are no longer cooling down
if ~isempty(tis.t1List)
    I = new_time - tis.t1Time < 10;
    % 				tis.t1CoolDown = tis.t1CoolDown - 1;
    % 				I = tis.t1CoolDown <= 0;
    tis.t1Time(I) = [];
    tis.t1List(:,I) = [];
end

% Find vertex pairs that are too close
D(~tis.connectivity) = NaN;
D(tis.parameters.fixed_verts,:) = NaN;
D = triu(D);
D( D==0 ) = NaN;
[I,J] = find( D < tis.parameters.t1Threshold );
for i = 1:numel(I)
    % Figure out which vertex-pairs are still in cooldown
    vID2transit = sort( vIDList([I(i) J(i)]) );
    match = bsxfun(@eq,vID2transit',[tis.t1List vID2transit']);
    match = logical(sum(match,1));
    % Perform T1 transitions
    if all(match <2)
        tis.t1Transition( tis.getVertices(vID2transit) );
    end
end

% Record energy
tis.energy = tis.get_energy;

end % move_vts