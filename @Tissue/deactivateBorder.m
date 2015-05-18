function deactivateBorder(tis,n)
% deactivateBorder:
% Deactivates the border (up to n order) cells
% USAGE: tis.deactivateBorder(n)
%
% @todo: implement orders n > 1

% By default only deactivate first layer
if nargin < 2, n = 1; end

if n>1, error('Haven''t implemented higher order :('); end

cellIDList = tis.cells.keys; cellIDList = [cellIDList{:}];
for ID = cellIDList
    if tis.numCellNeighbors( tis.cells(ID) ) < 6
        tis.deactivateCell( ID );
    end
end

end % deactivateBorder