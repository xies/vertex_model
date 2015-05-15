function varargout = solve_model(tis,ode_method,tspan,OUT_DIR,ode_opts)
%SOLVE_MODEL
% Uses the given ODE_METHOD to solve the model with curent
% Tissue as initial condition. And then assemble solver outputs
% using the assembler.
%
% Usage:
%   tis_init.solve_model( ...
%           ode_method,tspan,OUT_DIR,ode_opts);
%   tisArr = tis_init.solve_model( ...
%           ode_method,tspan,OUT_DIR,ode_opts);
%
% INPUT: ode_method - e.g. ODE23, ODE45
%        OUT_DIR - output directory for model (needs to be
%        clean)
%        tspan - time of solution
%        ode_opts (optional) - options for the ODE solver
%
% OUTPUT: optional: assembled tissue array
%
% See also: Tissue/evolve, Tissue/step, assemble_model

% Need to make sure everything is gone in OUT_DIR that might
% conflict.
s = what(OUT_DIR);
if ~isempty(s.mat)
    error(['Please clear the contents of ' OUT_DIR]);
end

csvwrite([OUT_DIR '/t1List.csv'], [NaN NaN NaN]);

verts = tis.vert_coords;
switch nargin
    case 4
        [T,Y] = ode_method(@(t,y) tis.step(t,y,OUT_DIR),tspan,verts);
    case 5
        [T,Y] = ode_method(@(t,y) tis.step(t,y,OUT_DIR),tspan,verts,ode_opts);
    otherwise
end
if nargout > 0
    varargout{1} = assemble_model(OUT_DIR);
end

csvwrite([OUT_DIR '/times.csv'],'T');
csvwrite([OUT_DIR '/vertices.csv'],'Y');

end % solve_model