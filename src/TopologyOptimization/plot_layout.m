function plot_layout(nel, v, varargin)
% Plots the layout of the topology optimization problem.
% Inputs:
%   nel: vector of size 2 x 1, the number of elements in each direction.
%   v: vector of size n x 1, the element densities or the signed distance
%       values of the level-set function.
%   mapFea2To: vector of size n x 1, the mapping from the FEA domain to the
%       TO domain (optional, default is []).

% Parse inputs
p = inputParser;
addOptional(p, 'mapFea2To', []);
parse(p, varargin{:});
mapFea2To = p.Results.mapFea2To;

% Pass from the FEA domain to the TO domain
if ~isempty(mapFea2To)
    v = v(mapFea2To);
end

nVals = numel(v);
nElements = prod(nel);

% Check the size of the input v
if nVals == nElements % the input v contains the element densities
    % Reshape the input v
    vMatrix = reshape(v, nel).';

    % Plot
    image(flipud(vMatrix), 'CDataMapping', 'scaled');
    colormap(gca, flip(gray));
    clim([0, 1])
    axis equal tight;
    box on;
    set(gca, 'FontSize', 20, 'xtick', [], 'ytick', []);
else % the input v contains the signed distance values of the level-set function
    % Reshape the input v
    vMatrix = reshape(v, nel + 1).';

    % Plot
    contourf(vMatrix, [0, 0], 'FaceColor', 'k');
    axis equal tight;
    box on;
    set(gca, 'XTick', [], 'YTick', [])
end
