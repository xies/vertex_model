function [vertlist] = get_vertices(bords)
% just gets the simple vertices by definition and returns the list

f = @(x) (x(2,2) && (sum(x(:)) >= 4));
lut = makelut(f, 3);

vertplot = bwlookup(bords, lut);

[y,x] = find(vertplot);

vertlist = [y(:),x(:)];

end