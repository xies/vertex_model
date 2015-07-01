function setContractilityModel(tis, contractions, cIDs)
% setContractilityModel - sets the value of contractility
% according to the given model.
%
% USAGE: tis.setContractilityModel( modelfun, params,cIDs)
%        tis.setContractilityModel( modelfun, params)
%                   - sets contractility in all ACTIVE cells
%
% INPUT: modelfuns - a function handle, should return the
%                   contractility if given the centroid
%                   locations (Nc x 2) and the parameters also
%                   passed to this function.
%
%        params - an array of parameters to pass to function
%
% Records the modelfun + params in tis.contraction_param

if nargin < 6
    cellsOI = tis.getActiveCells;
else
    cellsOI = tis.getCells( cIDs );
end

ct = cat(1,cellsOI.centroid);

% Only set spatial model at time = 0
if tis.t == 0
    
    % Synthesize modelfun and set the values
    C = synthesize_contractModels( ...
        ct,tis.t, contractions.spatial_model, contractions.spatial_params);
    tis.setContractility( C,[cellsOI.cellID] );
    tis.deactivateBorder;
    
    % Record modelfun + param
    tis.contractile_params = contractions;
    
    % Record initial contractility
    tis.contractile_params.initC = [tis.getActiveCells.contractility];
    
else % Use temporal model to evolve contractility
    
    C = synthesize_tempModels( tis, ...
        tis.contractile_params.temporal_model, ...
        tis.contractile_params.temporal_params );
    C = max(C,0);
    tis.setContractility( C,[cellsOI.cellID] );
    
end

end %setContractilityModel