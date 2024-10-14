function plot_backbone(omega, z, omegaTarget, zTarget)
% Plot the backbone curve in the physical domain (omega vs z).
% Additional points can be plotted by providing the omegaTarget and zTarget vectors.
% Inputs:
%   omega: vector, the natural frequencies [rad/s].
%   z: vector, the RMS values of the response [m].
%   omegaTarget: vector, the natural frequencies of the target points [rad/s]. Default is [].
%   zTarget: vector, the RMS values of the target points [m]. Default is [].

% Check inputs
if nargin == 2
    omegaTarget = [];
    zTarget = [];
end

% Plot backbone
hold on; grid on;
plot(omega, z, 'k', 'LineWidth', 2)

% Plot target points (if any)
if ~isempty(zTarget)
    plot(omegaTarget, zTarget, '.r', 'MarkerSize', 20)
    legend('Backbone', 'Target Points', 'Interpreter', 'latex')
end

% Decorate plot
xlabel('$\Omega$ [rad/s]', 'Interpreter','latex')
ylabel('$Z$ [m]', 'Interpreter','latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
axis square
box on
