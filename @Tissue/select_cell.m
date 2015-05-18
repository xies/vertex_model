function cellID = select_cell(tis,varargin)
%SELECT_CELL
% Shows the current tissue configuration and asks the user for a cell to
% select; outputs the cellID of the cell that contains that point(s).
% 
% USAGE: cellID = tis.select_cell;  (select single cell)
%        cellID = tis.select_cells(n);  (select n cells)
    
    if nargin > 1, n = varargin{1};
    else n = 1; end
    
    tis.draw;
    title('Select cell(s):')
    % Ask user to give query points
    [X,Y] = ginput(n);
    
    % Grab cell centroids
%     cx = tis.getCentroidX; cy = tis.getCentroidY;
%     cellIDList = tis.cells.keys; cellIDList = [cellIDList{:}];
    cellID = zeros(1,n);
    for j = 1:n
        x = X(j); y = Y(j);
        % Find the sq-distances from current point to all centroids
        sortedCells = tis.getCells.sortByDistance([x,y]);
        % Find the CellModel that contains point, from closest to farthest
        for c = sortedCells
            if c.contains( [y,x] ,tis)
                cellID(j) = c.cellID;
                break
            end
        end
    end

end