function setParameters(tis,varargin)
%setParameters - Sets the simluation/evolution parameters.
%
% USAGE: tissue =
%           tis.setParameters('parameter',value);
%   PARAMETERS:
%       'lineTension' - scalar or Nv x Nv matrix, no defaults
%       'lineAnisotropy' - anisotropic tension vertical v.
%                      horizontal junctions, default = 1
%       'areaElasticity' - scalar or Nc x 1 vector, no
%                       defaults
%       'perimeterElasticity' - scalar or Nc x 1 vector, no
%                       defaults
%       'connectivity' - string: see Matrix
%       'targetArea' - scalar or Nc x 1 vector, default is the
%                       mean value of current setting
%       'targetPerimeter' - default = 0
%       'dimensonless' - 1 for reduced parameter values,
%                       default = 0
%       'lengthScale' - for reducing parameters
%       'forceScale' - usually equal to linetension unless
%              there is special cases set or anisotropy
%       'viscosity' - scalar, no defaults
%       -- other parameters
%       'jitterSize' - for jittering vertex locations (px)
%                       default = 0
%       't1Threshold' -- default = 1 (px!)
%       'timeStep' - for Euler integration, default = 0.001
%       'frame_per_sec' - default = 1
%       'um_per_px' - default = 1
%
% xies@mit April 2015

p = myParameterParser;
Nb = tis.interfaces.length;
Nc = tis.cells.length;

% --- Dimensonality ---
% NO IDEA why matlab throws errors all over the place if you
% have a validator here...
addOptional(p,'dimensionless',false,@(x) isscalar(x)); %default = 0

% --- length scale ---
addOptional(p,'lengthScale',1,@isscalar); %default = 1

% --- force scale ---
addOptional(p,'forceScale',1,@isscalar); %default = 1

% --- Area elasticity ---
addRequired(p,'areaElasticity', ...
    @(x) isscalar(x) || (isvector(x) && numel(x) == Nc) );

% --- Line tension ---
addRequired(p,'lineTension',...
    @(x) isscalar(x) || (isvector(x) && numel(x) == Nb) );

% --- Line anisotropy ---
addOptional(p,'lineAnisotropy',1,@isscalar);

% --- Perimeter elasticity ---
addRequired(p,'perimElasticity', ...
    @(x) isscalar(x) || (isvector(x) && numel(x) == Nc) );

% --- Viscosity ---
addRequired(p,'viscosity',@isscalar);

% --- Connectivity ---
validConnectivityOpts = {'purse string'};
addRequired(p,'connectivity',...
    @(x) any(validatestring(x,validConnectivityOpts)));

% --- target area ---
c = tis.getCells; % default = mean of current areas
addOptional(p,'targetArea',mean([c.area]), ...
    @(x) isscalar(x) || (isvector(x) && numel(x) == Nc) );

% --- target perimeter ---
addOptional(p,'targetPerimeter',0, ... % default = 0
    @(x) isscalar(x) || (isvector(x) && numel(x) == Nc) );

% --- bookkeeping ---
addOptional(p,'jitterSize',0,@isscalar); %default = 0
addOptional(p,'stepSize',0.001,@isscalar); %default = 0.001
addOptional(p,'t1Threshold',2,@isscalar) % default = 1

% --- Units ---
addOptional(p,'um_per_px',1,@isscalar); %default = 1
addOptional(p,'dt_per_frame',1,@isscalar); %default = 1

parse(p,varargin{:});

tis.parameters = p.Results;
tis.parameters.fixed_verts = cellfun( ...
    @numel,{tis.getVertices.cellIDs}) < 3;

% Correct for um_per_px
tis.parameters.targetArea = ...
    tis.parameters.targetArea * tis.parameters.um_per_px^2;
tis.parameters.targetPerimeter = ...
    tis.parameters.targetPerimeter * tis.parameters.um_per_px;

if tis.parameters.dimensionless
    % Find dimensional parameters
    sigma = tis.parameters.forceScale; % force dimension
    lambda = tis.parameters.lengthScale; % length simension
    K_a = tis.parameters.areaElasticity * sigma / lambda^3;
    K_p = tis.parameters.perimElasticity * sigma / lambda;
    A0 = lambda^2 * 3 * sqrt(3) / 2;
else
    lambda = tis.parameters.lengthScale; % length simension
    sigma = tis.parameters.lineTension;
    K_a = tis.parameters.areaElasticity;
    K_p = tis.parameters.perimElasticity;
    A0 = tis.parameters.targetArea;
end

% Set Intarface-related parameters
tis.setLineTension;
if tis.parameters.lineAnisotropy ~= 1
    tis.setLineAnisotropy;
end

% Set CellModel parameters
tis.setAreaElasticity;
tis.setPerimElasticity;
tis.setTargetArea;

% Sanity check paramter regime
display('---Parameter check---')
display(['Normalized line tension = '...
    num2str( sigma / K_a / A0^3/2)])
display('(Should be ~0.12)')
display('-')
display(['Normalized perimeter tension = '...
    num2str( K_p / K_a / A0)])
display('(Should be ~0.04)')
display('---------------------');
display('')
display('---Force magnitud check---')
fMag = (- 6*sigma * lambda + K_a * A0^2 + K_p * (6*lambda)^2) ...
    * tis.cells.length;
display(['Force magnitude = ' num2str( fMag ) ])
display('-')
display(['Velocity magnitude = '...
    num2str( fMag * tis.parameters.stepSize / ...
    tis.parameters.viscosity / tis.parameters.um_per_px ) ])
display('---------------------');

% Update cells and interfaces
cellIDList = tis.cells.keys; cellIDList = [cellIDList{:}];
for i = 1:tis.cells.length
    tis.cells(cellIDList(i)) = ...
        tis.cells(cellIDList(i)).updateCell(tis);
end
bIDList = tis.interfaces.keys; bIDList = [bIDList{:}];
for i = 1:tis.interfaces.length
    tis.interfaces(bIDList(i)) = ...
        tis.interfaces(bIDList(i)).updateInterface(tis);
end

end % setParameters