function deactivateCell(tis, cellIDs)
% Deactivates cell(s)
% USAGE: tis.deactivateCell; (default = all cells)
%        tis.deactivateCell(IDs);
if nargin < 2 % If no cells specified, deactivate all cells
    cellIDs = tis.cells.keys;
    cellIDs = [cellIDs{:}];
end
for i = 1:numel(cellIDs)
    tis.cells( cellIDs(i) ) = ...
        tis.cells( cellIDs(i) ).deactivateCell;
end
end % deactivateCell