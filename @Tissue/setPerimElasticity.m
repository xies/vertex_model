function setPerimElasticity(tis,cIDs,alt_kp)
% Passes parameters from Tissue into CellModels
%
% Usage: setPerimElasticity(tis,cIDs,alt_kp)
%
% INPUT:
%   cIDs - by default does every cell
%   alt_kp - should match size of cIDs

if nargin < 2
    cIDs = tis.cells.keys;
    cIDs = [cIDs{:}];
    alt_kp = tis.parameters.perimElasticity;
end

num_cells = numel(cIDs);
% Update targetAreas
for i = 1:num_cells
    c = tis.cells( cIDs(i) );
    if isscalar(alt_kp)
        c.perimElasticity = alt_kp;
    else
        c.perimElasticity = alt_kp(i);
    end
    tis.cells( cIDs(i) ) = c;
end
end % setPerimElasticity