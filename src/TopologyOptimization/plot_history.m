function plot_history(history, legend_str)
% Plots the history of the optimization algorithm.
% Inputs:
%   history: matrix of size (n+1) x m, where n is the number of iterations
%            and m is the number of objectives + constraints.
%            The first column is the objective function and the rest are the
%            constraints.
%   legend_str: cell array of strings with the names of the objectives and
%               constraints (optional).

% Data
nPoints = find(~isnan(history(1, :)), 1, 'last');
iter = 1:nPoints;
plotData = history(:, 1:nPoints);

% Check inputs
nPlots = size(plotData, 1);
if nargin < 2
    legend_str = cell(1, nPlots);
    legend_str{1} = 'Objective Function';
    for ii = 2:nPlots
        legend_str{ii} = ['Constraint ', num2str(ii-1)];
    end
end

% Objective function
tiledlayout('flow', 'TileSpacing', 'compact');
nexttile;
hold on; grid on;
plot(iter, plotData(1, :), '-*k', 'LineWidth', 2);
xlabel('Iterations', 'Interpreter', 'latex');
title(legend_str{1}, 'Interpreter', 'latex');
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
axis tight
box on

% Constraints
for ii = 2:size(plotData, 1)
    nexttile;
    hold on; grid on;
    plot(iter, plotData(ii, :), '-*k', 'LineWidth', 2);
    xlabel('Iterations', 'Interpreter', 'latex');
    title(legend_str{ii}, 'Interpreter', 'latex');
    set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
    axis tight
    box on
end
