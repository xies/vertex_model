function C = time_of_start_linear_increase( tis, params )
%TIME_OF_START Starts nonzero contractility at specified time
%
% INPUT: PARAMS: 1x2 Cell array
%                params{1} - cellID
%                params{2} - time of start
%                params{3} - rate of change

p = zeros(1,numel(tis.getActiveCells));
p( params{1} ) = params{2};

C = max( tis.contractile_params.initC + ...
    params(1) * tis.t + params(2), ...
    0);

% C = ones(1,numel(C)) * max(C); % HACK!
C( tis.t < p ) = 0;

end