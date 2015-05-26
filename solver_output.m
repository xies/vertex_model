function status = solver_output(t,y,flag,varargin)
%PRINT
% Prints the current time to STDOUT
% For use with an odesolver.

% Status -- 0 is no problem / stop
status = 0;

display( ['Time = ' num2str(t)] );

end
