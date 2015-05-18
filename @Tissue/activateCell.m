function activateCell( tis, cellIDs, varargin)
%activateCell
% Set specified cells (IDs) to "active = 1". Used to keep track
% of who's "ventral" in the model.
% Usage: tis.activateCell(cellIDs)
%        tis.activateCell(cellIDs, alt_ka)
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
    nonActiveKA = varargin{1};
    inactiveCells = tis.getInactiveCells;
    num_cells = numel(inactiveCells);
    for i = 1:num_cells
        
        c = inactiveCells(i);
        c.areaElasticity = nonActiveKA;
        %                     edges = tis.getInterfaces( inactiveCells(i).bondIDs );
        tis.cells(c.cellID) = c;
        
    end
end

end % activateCell
