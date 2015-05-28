function setLineTension(tis,bIDList,L)
% Set the line tension of each bond in the tissue.
% Can specify each bond or have them be all the same
%
% USAGE: tis.setLineTension(bIDs,L)
%        tis.setLineTension - uses tis.parameters
%
% INPUT:
%   tis.parameters.lineTension:
%       scalar or Nv by Nv matrix with the same sparsity as
%       tis.connectivity

if nargin < 3
    bIDList = tis.interfaces.keys; bIDList = [bIDList{:}];
    L = tis.parameters.lineTension;
end

for i = 1:numel(bIDList)
    e = tis.interfaces( bIDList(i) );
    if isscalar(L)
        e.lineTension = L;
    else
        e.lineTension = L(i);
    end
    tis.interfaces( bIDList(i) ) = e;
end
if nargin > 1
    allIDs = tis.interfaces.keys; allIDs = [allIDs{:}];
    tis.parameters.lineTension( ismember(bIDList,allIDs) ) = L;
end

end % setLineTension