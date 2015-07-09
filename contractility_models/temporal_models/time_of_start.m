function C = time_of_start( tis, params )
%TIME_OF_START Starts nonzero contractility at specified time
%
% INPUT: PARAMS:

p = zeros(1,numel(tis.getActiveCells));
p( randi(numel(p),[1 floor(numel(p)/2)]) ) = params(1);

C = [tis.getActiveCells.contractility];
C = ones(1,numel(C)) * max(C); % HACK!
C( tis.t < p ) = 0;

end