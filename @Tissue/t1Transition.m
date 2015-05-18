function t1Transition( tis, vt )
%t1Transitions
%
% Make a T1 transition between given vertices. Not yet done!!!
%
% USAGE: tis.t1Transition( vt )
%
% INPUT: vt - two vertices to T1 transition

if numel(vt) ~= 2
    error('More than 2 vertices given!');
end
display(['T1 transition between ' num2str(vt(1).ID), ...
    ' and ' num2str(vt(2).ID) ])

% Swap cell ownership
% figure out which cells contain each vertex
cells1 = vt(1).cellIDs;
cells2 = vt(2).cellIDs;
% Figure out which bonds contain both vertices (should be only
% one!)
edge = tis.interfaces( intersect(vt(1).bondIDs, vt(2).bondIDs) );
other_edges = tis.getInterfaces( [vt.bondIDs] );

% figure out which cells contain both or just one
cells_both = tis.getCells( intersect(cells1,cells2) );
cells_single = ...
    tis.getCells( setxor(cells1,cells2) );

if numel(cells_both) > 2 || numel(cells_single) > 2
    display('T1 transition ill-defined!')
    keyboard % ill-defined!
end

% Rotate vertex positions to match the angle of cells_both
% and add both cells_single IDs to both vertices
midpt(1) = mean([vt.x]); midpt(2) = mean([vt.y]);
for i = 1:numel(vt)
    theta = cells_both.get_angle;
    vt(i) = vt(i).rotate(midpt,pi-theta);
    tis.vertices( vt(i).ID ) = vt(i);
    vt(i).cellIDs = union( vt(i).cellIDs, ...
        [cells_single.cellID]);
end

% Swap interface cell ownership and update tissue
edge.cIDs = [cells_single.cellID];
tis.interfaces(edge.ID) = edge.updateInterface(tis);

% Cells with only single vertex before now has both vertices
% and the bond
for i = 1:numel(cells_single)
    % Add vID to cell list
    cells_single(i).vIDs = union( cells_single(i).vIDs, [vt.ID]);
    % Add bondID to cell list
    cells_single(i).bondIDs = union( ...
        cells_single(i).bondIDs, edge.ID );
    % Put cell in Tissue
    tis.cells( cells_single(i).cellID ) = ...
        cells_single(i).updateCell(tis);
end

% Figure out which vertex is closer to each cells_both
D = cells_both(1).get_distance_to( [[vt.x]',[vt.y]'] );
[~,I] = max(D);
for i = 1:numel(cells_single)
    % Find the vertex that's farther away
    % remove the vID from cellModel
    ind = wrap(I+i-1,2);
    other_ind = wrap(I+i,2);
    cells_both(i).vIDs = setdiff(cells_both(i).vIDs, vt(ind).ID );
    % remove bondID from cellModel
    cells_both(i).bondIDs = setdiff(cells_both(i).bondIDs, edge.ID );
    
    % Get bonds that have this CellModel and remove the further
    % vertex and add the nearest vertex
    bIDs = intersect([other_edges.ID],cells_both(i).bondIDs);
    for j = 1:2
        this_bond = tis.interfaces( bIDs(j) );
        % Add/delete vID to INTERFACE
        this_bond.vIDs = setdiff( this_bond.vIDs, vt(ind).ID );
        this_bond.vIDs = union( this_bond.vIDs, vt(other_ind).ID );
        tis.interfaces( bIDs(j) ) = this_bond.updateInterface(tis);
        % Delete bondID from VERTEX and tis.CONNECTIVITY
        vt(ind).bondIDs = setdiff( vt(ind).bondIDs, this_bond.ID );
        % Add bondID from VERTEX
        vt(other_ind).bondIDs = unique( [vt(other_ind).bondIDs this_bond.ID] );
        
    end
    
    % remove cellID from vertex Model
    vt(ind).cellIDs = setdiff( vt(ind).cellIDs, cells_both(i).cellID );
    % Put vertex and cells back in tissue
    tis.vertices( vt(ind).ID ) = vt(ind);
    tis.vertices( vt(other_ind).ID ) = vt(other_ind);
    tis.cells( cells_both(i).cellID ) = ...
        cells_both(i).updateCell(tis);
end

%%% DEBUG!!!!!!!!!! %%%
%             if any(cellfun(@numel,{vt.bondIDs}) ~= 3)
%                 keyboard
%             end
%%%%%%%%

% Update matrices
tis.updateVertCoords;
tis.connect_interfaces('update');
tis.interVertDist = squareform(pdist( tis.vert_coords ));

% Keep track of T1 transitions
% Set cool down time to 10 steps
tis.t1List = cat(2,tis.t1List,sort([vt.ID])');
tis.t1Time = [tis.t1Time tis.t];



end % t1Transition