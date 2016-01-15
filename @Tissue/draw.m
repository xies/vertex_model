function I = draw(tis,varargin)
% Draws a single tissue in binary image. Can return just the
% outlines of cells, or also shade-in the active cells.
%
% USAGE: I = tis.draw;
%        I = tis.draw('showActive');
%           Highlights currentlly "active" cells
%        I = tis.draw('showContractile')
%           Highlights contractile cells and false-colors
%           contractility amount.
%        I = tis.draw('showVectors',V);
%           Shows vector field V over all vertices
%        I = tis.draw('showVectors',{v,ID});
%           Shows only single vector v at vertex == ID
%
% Cannot handle more than one tissue. Use tisArray.movie() for
% making a movie.

if numel(tis) > 1, error('Can only handle single tissue; use tis.movie() to show movie.'); end

%             um_per_px = tis.parameters.um_per_px;
%             x_axis = ((1:tis.Xs) - tis.Xs/2) * um_per_px;
%             y_axis = ((1:tis.Ys) - tis.Ys/2) * um_per_px;

I = zeros(tis.Xs,tis.Ys);
imagesc(I)
%             cellIDList = tis.cells.keys();
%             for i = 1:numel(cellIDList)
%                 I = I + tis.cells(cellIDList{i}).draw(tis);
%                 I = logical(I);
%             end
%             I = double(I)* 255;
hold on

% Highlight active cells
if any(strcmpi(varargin, 'showActive'))
    % Show active cells as filled-ins
    M = zeros(tis.Xs,tis.Ys);
    Acells = tis.getActiveCells;
    for i = 1:numel(Acells)
        M = M + Acells(i).drawMask(tis);
    end
    M = M * 50;
    I = I + M;
end

% Highlight active cells
if any(strcmpi(varargin, 'showContractile'))
    % Show active cells as filled-ins, with value proportional
    % to contractility
    M = zeros(tis.Xs,tis.Ys);
    Acells = tis.getActiveCells;
    C = [Acells.contractility];
    cmax = max(C);
    %                 cmax = 0.05;
    for i = 1:numel(Acells)
        M = M + Acells(i).drawMask(tis) * ...
            Acells(i).contractility * 150 / cmax;
    end
    I = I + M;
end

% Highlight specific cells
ind = find( strcmpi(varargin,'showCellID') );
if ~isempty(ind)
    M = zeros(tis.Xs,tis.Ys);
    IDs = varargin{ind+1};
    Acells = tis.getCells(IDs);
    for i = 1:numel(Acells)
        M = M + Acells(i).drawMask(tis);
    end
    I = I + M * 50;
end

hold off, imagesc(I), axis equal; colormap hot, colorbar hot, caxis([0 255])
%             hold off, imagesc(y_axis,x_axis,I), axis equal;
hold on, tis.getInterfaces.draw(tis);

% Highlight vertices
if any(strcmpi(varargin, 'showVertices'))
    hold on
    vcoords = tis.vert_coords;
    scatter(vcoords(:,2),vcoords(:,1),'w');
end

% Display vector field
ind = find( strcmpi(varargin,'showVectors') );
if ~isempty(ind)
    if numel(varargin) < ind + 1, error('Need vector to draw'); end
    V = varargin{ind+1};
    % If input is not a cell object, then display all vectors
    if ~iscell(V)
        figure(1)
        hold on;
        quiver(...
            tis.vert_coords(:,2), ...
            tis.vert_coords(:,1), ...
            V(:,2), ...
            V(:,1), ...
            0,'r-');
    else
        % If it's a cell obj, then only the given vertex will
        % have a vector over it
        v = V{1};
        ID = V{2};
        if numel(ID) ~= size(v,1);
            error('# of vectors should equal # of origins')
        end
        hold on
        quiver(tis.vert_coords(ID,2),tis.vert_coords(ID,1), ...
            v(:,2),v(:,1),0,'w-');
    end
end

axis off
hold off

end % DRAW