function flag = isValid(tis)
% Run self-consistency tests on current tissue state
% Checks the following:
%   1) vertices and vert_coords match and are ordered correctly
%   2) vertices and the unique set of cell-owned vertices also
%      match (not ordered)
%   3) dimensions of connectivity and interVertDist are
%      consistent with vert_coords

flag = 1;
% 1) vertices and vert_coords match
vx = [tis.getVertices.x]; vy = [tis.getVertices.y];
flag = flag & all(all(cat(2,vx',vy') == tis.vert_coords));

% 2a) Each cell's list of vertices match vertex cellIDs
vIDs = {tis.getCells.vIDs};
v = cellfun(@tis.getVertices,vIDs,'UniformOutput',0);
v = cellfun(@(x) cat(1,x.cellIDs),v,'uniformoutput',0);
flag = flag & all(cellfun(@(x,y) any(x == y), v, {tis.getCells.cellID}));

% 2b) Each vertex's list of cells match cell vIDs
cellIDs = {tis.getVertices.cellIDs};
c = cellfun(@tis.getCells,cellIDs,'UniformOutput',0);
c = cellfun(@(x) cat(2,x.vIDs),c,'uniformoutput',0);
flag = flag & all(cellfun(@(x,y) any(x == y), c, {tis.getVertices.ID}));

% 3a) Each cell's list of edges match Interface cellIDs
bIDs = {tis.getCells.bondIDs};
e = cellfun(@tis.getInterfaces,bIDs,'UniformOutput',0);
e = cellfun(@(x) cat(2,x.cIDs),e,'uniformoutput',0);
flag = flag & all(cellfun(@(x,y) any(x == y), e, {tis.getCells.cellID}));

% 3b) Each bond's list of cells match cell bIDs
cellIDs = {tis.getInterfaces.cIDs};
c = cellfun(@tis.getCells,cellIDs,'UniformOutput',0);
c = cellfun(@(x) cat(2,x.bondIDs),c,'uniformoutput',0);
flag = flag & all(cellfun(@(x,y) any(x == y), c, {tis.getInterfaces.ID}));

% 4a) Each bond's list of vertices match vertex bIDs
vIDs = {tis.getInterfaces.vIDs};
v = cellfun(@tis.getVertices,vIDs,'UniformOutput',0);
v = cellfun(@(x) cat(2,x.bondIDs),v,'uniformoutput',0);
flag = flag & all(cellfun(@(x,y) any(x == y), v, {tis.getInterfaces.ID}));

% 4b) Each vertex's list of bonds match bond vIDs
bIDs = {tis.getVertices.bondIDs};
e = cellfun(@tis.getInterfaces,bIDs,'UniformOutput',0);
e = cellfun(@(x) cat(2,x.vIDs),e,'uniformoutput',0);
flag = flag & all(cellfun(@(x,y) any(x == y), e, {tis.getVertices.ID}));

end % isValid