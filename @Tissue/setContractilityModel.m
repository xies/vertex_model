function setContractilityModel(tis,modelfun,params,cIDs)
% setContractilityModel - sets the value of contractility
% according to the given model.
%
% USAGE: tis.setContractilityModel( modelfun, params,cIDs)
%        tis.setContractilityModel( modelfun, params)
%                   - sets contractility in all ACTIVE cells
%
% INPUT: modelfun - a function handle, should return the
%                   contractility if given the centroid
%                   locations (Nc x 2) and the parameters also
%                   passed to this function.
%
%        params - an array of parameters to pass to function
%
% Records the modelfun + params in tis.contraction_param

if nargin < 4
    cellsOI = tis.getActiveCells;
else
    cellsOI = tis.getCells( cIDs );
end

% Evaluate modelfun and set the values
ct = cat(1,cellsOI.centroid);
C = modelfun(ct,tis.t,params);
tis.setContractility( C,[cellsOI.cellID] );
tis.deactivateBorder;

% Record modelfun + param
tis.contractile_params.model = modelfun;
tis.contractile_params.params = params;

end %setContractilityModel