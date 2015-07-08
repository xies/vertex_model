function C = time_of_start( tis, params )
%TIME_OF_START Starts nonzero contractility at specified time
%
% INPUT: PARAMS:

C = [tis.getActiveCells.contractility];
C = ones(1,numel(C)) * max(C); % HACK!
C( tis.t < params ) = 0;

end