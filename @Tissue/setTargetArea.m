function setTargetArea(tis,cIDs,alt_area)
% Passes parameters from Tissue into CellModels
%
% Usage: setTargetArea(tis,cIDs,alt_area)
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
end % setTargetArea