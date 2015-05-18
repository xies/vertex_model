function setLineAnisotropy(tis,bIDList,a)
% Set the line tension anisotropically according to
% orientation.
%
% Anisotropy value is the same across specified cells
%
% USAGE: tis.setLineAnisotropy(bIDs,a)
%        tis.setLineAnisotropy - use tis.parameters as
%           default value

sigma = tis.parameters.forceScale;
if nargin < 3
    a = tis.parameters.lineAnisotropy;
    bIDList = tis.interfaces.keys; bIDList = [bIDList{:}];
end

if a == 1, return; end % Do nothing if a is unit

new_tensions = zeros(1,numel(bIDList));
for i = 1:numel(bIDList)
    e = tis.interfaces( bIDList(i) );
    theta = e.angle;
    % Scale horizontal junctions by factor
    if theta < pi/4 || theta > 3*pi/4
        e.tension = sigma*a;
        % set new junctions
        tis.interfaces( bIDList(i) ) = e;
    end
    new_tensions(i) = e.tension;
end

% Update parameter list
allIDs = tis.interfaces.keys; allIDs = [allIDs{:}];
if isscalar( tis.parameters.lineTension )
    tis.parameters.lineTension = ...
        tis.parameters.lineTension*ones(1,numel(allIDs));
end
tis.parameters.lineTension( ismember(bIDList,allIDs) ) = new_tensions;

end % setLineAnisotropy