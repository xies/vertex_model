function C = time_of_start( tis, params )
%TIME_OF_START Starts nonzero contractility at specified time
%
% INPUT: PARAMS: 1x2 Cell array
%        params{1} - cellID

p = zeros(1,numel(tis.getActiveCells));
p( params{1} ) = params{2};

C = [tis.getActiveCells.contractility];
C = ones(1,numel(C)) * max(C); % HACK!
C( tis.t < p ) = 0;

end