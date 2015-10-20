function C = linear_increase( tis, params )
%LINEAR_INCREASE

C = max( tis.contractile_params.initC + ...
    params(1) * tis.t + params(2), ...
    0);
%     tis.parameters.areaElasticity * 0.1 );

end