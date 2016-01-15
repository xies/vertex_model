function [centroids,vertices,regions] = create_random_grid(Ncells,Nsmooth)
% Use Lloyd's method to generate random polygonal grids with variable
% smoothness (higher Nsmooth -> more like regular hexagon)
%
% USAGE: [centroids,vertices,regions] =
%           create_random_grid(Ncells,Nsmooth)
%
% xies@mit Jan 2016

maxN = 2^8;
boundary = [0 0; 0 maxN; maxN maxN; maxN 0];

centroids = rand(Ncells,2) * maxN;
% plot(centroids(:,1),centroids(:,2),'ro');

for t = 1:Nsmooth % The higher Nsmooth is, the more regular the grid is
    
    % Find Voronoi tesselation by delaunay triangulation
%     DT = delaunayTriangulation(centroids);
    [V,R] = voronoin( centroids );
    
    % Find regions on the convex hull
    ind2trim = find_on_edge(V,R,maxN);
    
    % Find centroids of tesselations unless on chull
    for i = setdiff(1:Ncells,ind2trim)
        centroids(i,:) = mean( V(R{i},:) );
    end
    
    for i = ind2trim
        I = R{i}; I( I == 1 ) = [];
        vts = intersectionHull('vert',boundary,'vert',V(I,:));
        centroids(i,:) = mean( vts.vert );
    end
    
end

% hold on
% plot(centroids(:,1),centroids(:,2),'bo');

% Give final vertices
[V,R] = voronoin(centroids);
vertices = V;

% Generate labeled binary image with regions
regions = false(maxN);
for i = 1:numel(R)
    I = R{i}; I(I==1) = [];
    regions = regions | bwmorph(poly2mask(V(I,1),V(I,2),2^8,2^8),'thin',1);
end
regions = bwmorph(regions,'thicken',Inf);
regions = bwlabel(regions);

% get rid of out of bound vertices
vertices(1,:) = [];
vertices(any(vertices > maxN,2),:) = [];
vertices(any(vertices < 0,2),:) = [];

end

function ind2trim = find_on_edge(V,R,maxN)
I = zeros(1,numel(R));
for i = 1:numel(R)
    I(i) = any(any(V(R{i},:) < 0 | V(R{i},:) > maxN));
end
ind2trim = find(I);

end