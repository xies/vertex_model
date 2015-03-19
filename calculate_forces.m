function [bigvert biggrad stress] = calculate_forces(verts, p)

% do this by cutting each strand and looking at the gradient
% a "cut" is equivalent to knocking out a 1 in the connectivity matrix
% this will give the forces for all the strands, and we can draw
% a quiver plot
% and use the formula for maximum strain (represent with color???)

preferred_distances = p.preferred_distances;
preferred_distances(~p.connectivity) = NaN;
preferred_distances(logical(eye(size(preferred_distances)))) = NaN; 

spring_constants = p.spring_constants;

stress = zeros(p.numV, 1);
bigvert = [];
biggrad = [];
distances = compute_distances(verts, p.numV);
initial_grad = compute_gradient(distances, preferred_distances, spring_constants, verts);
for f = 1:p.numV
    myrow = preferred_distances(f, :);
    % all the vertices that this vertex is connected to
    connects = find(~isnan(myrow));
    nconn = length(connects);

    tempgrad = [];
    for arrow = 1:nconn
        preferred_distances_f = preferred_distances;
        % the connection that will be knocked out
        knockout = connects(arrow);
        preferred_distances_f(f, knockout) = NaN;
        preferred_distances_f(knockout, f) = NaN;
        
        grad = compute_gradient_one_element(f, distances, preferred_distances_f, spring_constants, verts);
%             test = compute_gradient(distances, preferred_distances, verts);

        force = --grad;
        % when connection is cut, force points away
        % we want it to point towards, so another - sign

        % calculate what the force would be without a knockout- this will be 0 for
        % all except the fixed verts. then subtracts this from force to get the contribution
        % just from that pull

        initial_force = --initial_grad(f, :);

        net_force = force - initial_force;
        tempgrad = [tempgrad; net_force];


        % put the vertex coordinate shifted a bit
        % towards where the connection is
        % so that you can see which arrow is which
        sf = 100;  % the fraction to shift by
        placecoord = 1/sf*((sf-1)*verts(f, :) + verts(knockout, :));
        bigvert = [bigvert; placecoord];
    end

    TOL = 0.01;
    % check that the forces add to zero
    % but only for non-fixed cells (the first "if")
    if sum(abs(initial_force) > TOL) == 0     
        if abs(sum(tempgrad(:))) > TOL
            disp('Error! Forces do not add to zero');
            initial_force
            tempgrad
            keyboard
        end
    end


    % add entries to bigvert, to go into quiver
%     bigvert = [bigvert; repmat(verts(f, :), nconn, 1)];
    biggrad = [biggrad; tempgrad];

    stress(f) = tension(tempgrad);
end