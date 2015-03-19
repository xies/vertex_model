% rounds the elements of "in" to the nearest "n"
% round_n(in, 1) = round(in)
function out = round_n(in, n)

out = round(in / n) * n;


