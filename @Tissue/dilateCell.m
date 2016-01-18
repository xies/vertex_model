function dilateCell(tis,cIDs,scale)
%dilateCell Isotropically transform given cells (dilate or shrink)
% If multiple cellIDs are given, will do it in array order.
% 
% USAGE: tis.dilateCell( cIDs, scale );
%
% xies@mit.edu Jan 2016

num_cells = numel(cIDs);

for i = 1:num_cells
    cell = tis.getCells( cIDs );
    vts = tis.getVertices( cell.vIDs );
    vx = [vts.x];
    vy = [vts.y];
    
    new_vx = (vx - cell.centroid(1)) * scale + cell.centroid(1);
    new_vy = (vy - cell.centroid(2)) * scale + cell.centroid(2);
    
    for j = 1:numel(cell.vIDs)
        tis.moveVertex( cell.vIDs(j), [new_vx(j), new_vy(j)] );
    end
    
end

end