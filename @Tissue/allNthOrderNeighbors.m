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