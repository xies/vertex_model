function E = get_energy(tis)
% GET_ENERGY Returns the current energy of the system by adding
% the energy of all cells and interfaces.
%
% USAGE: E = get_energy(tis)

% Interface energy (line tension)
e = tis.interfaces.values; e = [e{:}];
bondTerm = sum([e.energy]);
% Cell energy (area and perim elasticity)
c = tis.cells.values(); c = [c{:}];
cellTerm = sum([c.energy]);

E = bondTerm + cellTerm;
end % get_energy