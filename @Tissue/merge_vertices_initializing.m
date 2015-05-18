function vcoords = merge_vertices_initializing(~,vcoords, ...
    merge_threshold_in_px)
% Merge vertices that are closer than the specified threshold
%
% USAGE:
% vcoords = tis.merge_vertices_initializing( vcoords,threshold)
%
% Use only from constructor

num_vertices = size(vcoords);
% Merge vertices which are super close to each other
vertDist = squareform(pdist(vcoords));
vertDist( logical(eye(num_vertices)) ) = NaN;

% Loop merge function until no pairs of vertices are closer than
% the threshold
while any(any(vertDist <= merge_threshold_in_px))
    
    % Find a set of vertices to merge
    [I,J] = find(vertDist <= merge_threshold_in_px,1,'first');
    tobeMergedInd = [I,find(vertDist(I,:) <= merge_threshold_in_px)];
    mergedV = mean( vcoords(tobeMergedInd,:) );
    
    % Delete old ummerged vertices
    vcoords(tobeMergedInd,:) = [];
    % Add new merged vertex
    vcoords = cat(1,vcoords,mergedV );
    
    % update the distance maps
    num_vertices = size(vcoords,1);
    vertDist = squareform(pdist(vcoords));
    vertDist( logical(eye(num_vertices)) ) = NaN;
    
end
end % merge_vertices_initializing