classdef Interface
    %INTERFACE
    %
    % Properties:
    %   ID - key to Map held by Tissue
    %   lineTension - energy per length (surface tension)
    %   vIDs - IDs of 2 Vertex objects connected to this Interface
    %   cIDs - the 2 CellModel objects touching this Interface
    %   angle - angle of interface w.r.t x-axis
    %   energy - sigma * length energy contribution of interface
    %   length - length in microns
    %
    % Methods:
    % --- Constructor ---
    %   Interface - constructor
    %
    % --- Measurements ---
    %	get_length
    %   get_angle
    %   get_energy
    %   updateInterface - update interface given new Tissue object
    %
    % --- Visualize ---
    %   draw - draw interface arrays as matlab LINE objects
    %   drawColor - draw colored by energy
    %
    % SEE ALSO: Tissue, CellModel, Vertex
    
    properties
        ID
        lineTension
        vIDs
        cIDs
        angle
        energy
        length
    end
    methods
        
        % --- Constructor ---
        function obj = Interface(ID,vts,cellIDs,tis)
            if nargin > 0
                if numel(vts) ~= 2, error('Requires 2 vertices as input'); end
                if vts(1) == vts(2), error('Requires 2 vertices as input'); end
                obj.ID = ID;
                obj.vIDs = [vts.ID];
                obj.cIDs = cellIDs;
                
                % --- Set measurements ---
                obj.angle = obj.get_angle(tis);
                obj.length = obj.get_length(tis);
                obj.lineTension = 1;
            end
        end
        
        % --- Measurements ---
        function l = get_length(bond,tis)
            % Gets the length of current edge
            vts = tis.getVertices( bond.vIDs );
            if ~isempty(tis.parameters)
                x = tis.parameters.um_per_px^2;
            else
                x = 1;
            end
            l = vts(1).distance(vts(2)) * x;
        end
        function theta = get_angle(bond,tis)
            % Returns the angle to x-axis in the I + II quadrants (0 to pi)
            vts = tis.getVertices( bond.vIDs );
            theta = atan2(diff([vts.x]),diff([vts.y] ));
        end
        function E = get_energy(bond,tis)
            % Returns the line tension energy value at this interface
            E = bond.length * bond.lineTension / ...
                tis.parameters.forceScale / tis.parameters.lengthScale;
        end
        
        function bond = updateInterface(bond,tis)
            % Update the interface length and energy according to the
            % current tissue configuration
            %
            % USAGE: bond = bond.updateInterface(tis)
            %
            bond.length = bond.get_length(tis);
            bond.energy = bond.get_energy(tis);
        end
        
        %  --- Visualize ---
        function draw(bonds,tis)
            for i = 1:numel(bonds)
                v = tis.getVertices(bonds(i).vIDs);
                v.line();
%                 v.draw()
                hold on
            end
            hold off
        end
        
        function drawColor(bonds,tis)
            tmax = max([bonds.energy]);
            trange = linspace(0,tmax,64);
            C = jet;
            for i = 1:numel(bonds)
                v = tis.getVertices(bonds(i).vIDs);
                v.line( C(findnearest(bonds(i).energy,trange),:) );
                hold on
            end
            hold off
        end
    end
    
end