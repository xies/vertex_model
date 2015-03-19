function draw_connectivity_matrix_fast(connect, list, col)
%connect is a connectivity matrix
%list is a list of centroids or vertices
%col is the color

[ind1 ind2] = find(triu(connect));

for i = 1:length(ind1)      
    line([list(ind1(i),2) list(ind2(i),2)], [list(ind1(i),1) list(ind2(i),1)], 'Color', col);
end
