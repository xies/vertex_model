function [coid regions] = get_cents(bords)
% This function gets the centroid position of the hexagon (center of mass of the 
% hexagon as well as the regions of the hexagons
regions = bwlabel(logical(1-bords), 4);
cent = regionprops(regions,'Centroid');
coid = [cent.Centroid];
coid = reshape(coid,[2 length(coid)/2]).';
coid = round(coid);

% flip x and y
coid = fliplr(coid);

%gets rid of the background
coid = removerows(coid, 1); 
%fix region array
regions = regions - 1;

end