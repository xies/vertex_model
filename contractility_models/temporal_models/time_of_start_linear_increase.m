function C = time_of_start_linear_increase( tis, params )
%TIME_OF_START Starts nonzero contractility at specified time
%
% INPUT: PARAMS: 1x2 Cell array
%                params{1} - cellID
%                params{2} - time of start
%                params{3} - rate of change

p = zeros(1,numel(tis.getActiveCells));
p( params{1} ) = params{2};

C = tis.contractile_params.initC + ...
    params{3} * tis.t;

C( p ) = C( p ) - params{2}*params{3};
C = max(C,0);

end