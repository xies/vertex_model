classdef Interface
    
    properties
        ID
        tension
        vIDs
        cIDs
        angle
    end
    methods
        
        % --- Constructor ---
        function obj = Interface(ID,vts,cellIDs)
            if nargin > 0
                if numel(vts) ~= 2, error('Requires 2 vertices as input'); end
                if vts(1) == vts(2), error('Requires 2 vertices as input'); end
                obj.ID = ID;
                obj.vIDs = [vts.ID];
                obj.cIDs = cellIDs;
%                 obj.coords = cat(2,[vts.x]',[vts.y]');
                
                theta = atan2(diff([vts.y]),diff([vts.x]));
                obj.angle = abs(theta);
            end
        end
        
        
        %  --- Visualize ---
        function draw(bonds,tis)
            tmax = max([bonds.tension]);
            trange = linspace(0,tmax,64);
            C = hot;
            for i = 1:numel(bonds)
                v = tis.getVertices(bonds(i).vIDs);
                v.line( C(findnearest(bonds(i).tension,trange),:) );
                hold on
            end
            hold off
        end
    end
    
end