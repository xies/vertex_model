function C = synthesize_tempModels(tis,modelFuns,modelParams)
%SYNTHESIZE_TEMPMODELS
% Updates the currently "active" cells according to the temporal update model
% given.
%
% C = synthesize_tempModels(tis,modelFuns,modelParams)

numFcn = numel(modelFuns);
numCells = numel(tis.getActiveCells);
C = zeros(numCells,1);
for i = 1:numFcn
    p{1} = tis; p{2} = modelParams{i};
    C = C + ensure_column( feval( modelFuns{i}, p{:} ) );
end

end