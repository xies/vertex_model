function moveVertex(tis,vID,new_coord)
% Move a single vertex
tis.vertices(vID) = tis.vertices(vID).move(new_coord);
tis.updateVertCoords;

v = tis.vertices(vID);

% Update cells associated with this vertex
cContainV = v.cellIDs;
for c = cContainV'
    tis.cells( c ) = ...
        tis.cells(c).updateCell( tis );
end

% Update edges associated with this vertex
eContainsV = v.bondIDs;
for e = ensure_row(eContainsV)
    tis.interfaces( e) = ...
        tis.interfaces(e).updateInterface(tis);
end

end