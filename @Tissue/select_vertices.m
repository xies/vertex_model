function vID = select_vertices(tis,varargin)
%SELECT_VERTICES
% Shows the current tissue configuration and asks the user for a query
% point(s) and returns vID of vertex closest to point(s).
% 
% USAGE: vID = tis.select_vertices;  (select single vertex)
%        vID = tis.select_vertices(n);  (select n vertex)
    
    if nargin > 1, n = varargin{1};
    else n = 1; end
    
    tis.draw;
    title('Select cell(s):')
    % Ask user to give query points
    [X,Y] = ginput(n);
    
    % Grab cell centroids
%     cx = tis.getCentroidX; cy = tis.getCentroidY;
%     cellIDList = tis.cells.keys; cellIDList = [cellIDList{:}];
    vID = zeros(1,n);
    for j = 1:n
        x = X(j); y = Y(j);
        % Find the sq-distances from current point to all centroids
        sortedVts = tis.getVertices.sortByDistance([x,y]);
        % Find the CellModel that contains point, from closest to farthest
        vID(j) = sortedVts(1).ID;
    end

end