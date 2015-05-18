function setAreaElasticity(tis,cIDs,alt_ka)
% Passes parameters from Tissue into CellModels
%
% Usage: setAreaElasticity(tis,cIDs,alt_ka)
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