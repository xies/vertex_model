function setContractility(tis,C,cIDs)
%setContractility
% Directly sets the active contractility coefficient
% Requires all cells to have its own specified contractility
%
% USAGE: tis.setContractility( C , cIDs );
%        tis.setContractility( C , cIDs );
% INPUT: tis - tissue
%        C - vector of contractility values
%        cIDs - by default sets all cells

if nargin < 3
    cIDs = tis.cells.keys; cIDs = [cIDs{:}];
end

if numel(cIDs) ~= numel( C )
    error('Number of contractility coeff and number of cells don''t match');
end

for i = 1:numel(cIDs)
    if tis.cells(cIDs(i)).isActive
        % set cell contractility iff it's an active cell
        tis.cells( cIDs(i) ) = ...
            tis.cells( cIDs(i) ).setContractility(C(i));
    end
end
end % setContractility