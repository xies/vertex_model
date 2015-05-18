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