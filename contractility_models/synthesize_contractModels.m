function C = synthesize_contractModels(ct,t,modelFuns,modelParams)
% SYNTHESIZE_CONTRACTMODELS
% Returns the sum of several different contractility models (e.g. spatial
% gradient + noise), if given the model function handles and parameters.
%
% USAGE: C = synthesize_contractModels(c,t,modelFuns,modelParams)
%
% INPUT: ct - current cell centroids (size: Nc x 2)
%        t - current simulation time
%        modelFuns - a cell array of model handles
%        modelParams - cell array of model parameters, order corresponds to
%           modelFuns
% OUTPUT: C - contractility array (same size as ct)

numFcn = numel(modelFuns);
C = zeros(1,size(ct,1));
for i = 1:numFcn
    p{1} = ct; p{2} = t; p{3} = modelParams{i};
    C = C + feval( modelFuns{i}, p{:} );
end

end