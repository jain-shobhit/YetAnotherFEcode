function plot_layout_3d(dim, coord, d, varargin)
    % Plots the 3d layout of the topology optimization problem.
    % Inputs:
    %   dim: the dimensions of the design domain.
    %   coord: the n x 3 matrix of coordinates.
    %   d: the density vector of length n.
    %   dLim: the threshold for the density (optional, default 0.5).

    % Parse the inputs
    p = inputParser;
    addOptional(p, "dLim", 0.5);
    parse(p, varargin{:});
    dLim = p.Results.dLim;

    % Extract the coordinates of the elements with density above the threshold
    coordFilt = coord(d > dLim, :);
    scatter3(coordFilt(:, 1), coordFilt(:, 2), coordFilt(:, 3), 50, d(d > dLim), "filled");

    % Add colormap
    colormap(gca, flip(gray));
    clim([0, 1])

    % Decorate the plot
    hold on; axis equal; box on;
    xlabel("X", "Interpreter", "latex")
    ylabel("Y", "Interpreter", "latex")
    zlabel("Z", "Interpreter", "latex")
    xlim([0.0, dim(1)])
    ylim([0.0, dim(2)])
    zlim([0.0, dim(3)])
    set(gca, "FontSize", 20, "XTick", [], "YTick", [], "ZTick", [])
end
