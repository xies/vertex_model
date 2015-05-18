function dy = step(tis,t,y,OUT_DIR)
%STEP
% Wrapper for move_vts function for use with native ODE solver
% like ODE23 or ODE45. Call using an anonymous function to pass
% the tissue and the directory settings. Will save each
% iteration as a .mat file if given a directory.
%
% USAGE: dy = tis.step(t,y,OUT_DIR)
%
% INPUT: tis - Tissue
%        t - current time
%        y - current vertices
%        OUT_DIR (optional) - output directory
%
% SEE ALSO: ODE23, ASSEMBLE_MODEL, RUN_MODEL

y = reshape(y(:),[numel(y)/2 2]);
tis.move_vts(y,t);
dy = tis.get_force / tis.parameters.viscosity;
dy = dy(:);

if nargin > 3
    save([OUT_DIR '/model_t_' num2str(t) '.mat'],'tis');
end

end % step