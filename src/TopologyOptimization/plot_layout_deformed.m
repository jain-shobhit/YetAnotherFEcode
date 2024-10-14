function plot_layout_deformed(nel, d, nodes, u, mapFea2To)
% Plots the layout of the topology optimization problem.
% Inputs:
%   nel: vector of size 2 x 1, the number of elements in each direction.
%   d: vector of size n x 1, the element densities.
%   nodes: matrix of size n x 2, the coordinates of the nodes.
%   u: vector of size 2 * n x 1, the displacement field.
%   mapFea2To: vector of size n x 1, the mapping from the FEA domain to the TO domain (optional).

% Check inputs
if nargin < 5
    mapFea2To = [];
end

% Pass from the FEA domain to the TO domain
if ~isempty(mapFea2To)
    d = d(mapFea2To);
end

% Reshape the input d
dMatrix = reshape(d, nel).';

% Scale the displacement
uLim = 0.5 * max(nel);
uMax = max(abs(u));
u = u / uMax * uLim;

% Extract displacement
uX = u(1:2:end);
uY = u(2:2:end);

% Sort nodes along X
[~, idxSort] = sort(nodes(:, 1));
nodes = nodes(idxSort, :);
uX = uX(idxSort);
uY = uY(idxSort);

% Sort nodes along Y
[~, idxSort] = sort(nodes(:, 2));
nodes = nodes(idxSort, :);
uX = uX(idxSort);
uY = uY(idxSort);

% Extract node coordinates
nodesX = reshape(nodes(:, 1) + uX, nel + 1).';
nodesY = reshape(nodes(:, 2) + uY, nel + 1).';

% Plot
surf(nodesX, nodesY, ones(size(nodesX)), dMatrix, 'LineStyle', 'none');
view(2);
colormap(gca, flip(gray));
clim([0, 1])
axis equal tight off;
set(gca, 'FontSize', 20, 'xtick', [], 'ytick', []);
