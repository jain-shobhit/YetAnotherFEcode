clear; clc; close all;
addpath(genpath(fullfile('..', 'TO-LSM-functions')))

load("beam_linear_pUpdate.mat")
% load("beam_hardening_pUpdate.mat")
% load("beam_softening_pUpdate.mat")

% History
figure();
plot_history(history);

% Initial layout
figure();
plot_layout(to.nel, densHistory(:, 1), to.mapFea2To);
% title('Initial Layout', 'Interpreter', 'latex');

% Optimal layout
figure();
plot_layout(to.nel, to.d_proj, to.mapFea2To);
title('Optimal Layout', 'Interpreter', 'latex');

% Evaluate gamma
T2 = myAssembly.tensor_uniform('T2', [nDOFs, nDOFs, nDOFs], [2, 3], 'weights', to.d_simp);
T3 = myAssembly.tensor_uniform('T3', [nDOFs, nDOFs, nDOFs, nDOFs], [2, 3, 4], 'weights', to.d_simp);
T2c = myAssembly.constrain_tensor(T2);
T3c = myAssembly.constrain_tensor(T3);
g = gamma_value(Mc, Kc, T2c, T3c, omega, uc);
fprintf('\nFinal gamma: %.4e\n\n', g)

% Get the target dof
targetNode = find_node(lx, ly/2, [], nodes);
targetDofs = get_index(targetNode, myMesh.nDOFPerNode);
targetDof = myAssembly.free2constrained_index(targetDofs(2));

% Compute rho from the target physical displacement
wTarget = 10;
rhoTarget = compute_rho(wTarget, Mc, Kc, T2c, T3c, omega, uc, targetDof);
% rhoTarget = 1e3;

% Compute manifold and backbone
rho = linspace(0, rhoTarget, 51);
[Omega, theta, wManifold, wBB] = compute_manifold(Mc, Kc, T2c, T3c, omega, uc, rho, targetDof);

% Plot manifold
figure
hold on; grid on; box on;
surf(rho.' * cos(theta), rho.' * sin(theta), wManifold, ...
    'FaceColor', [0.0745 0.6235 1], 'EdgeColor', 'none', 'FaceAlpha', 0.7)
for ii = [20, 40]
    plot3(rho(ii) * cos(theta), rho(ii) * sin(theta), wManifold(ii, :), ...
        'r.-', 'LineWidth', 2)
end
L = light;
L.Position = [0 0 1];
lighting gouraud
xlabel('$\rho \cos(\theta)$', 'Interpreter', 'latex')
ylabel('$\rho \sin(\theta)$', 'Interpreter', 'latex')
zlabel('$w$ [$\mathrm{\mu}$m]', 'Interpreter', 'latex')
title('Manifold', 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
axis square tight
box on
view(45,15)

% Plot backbone in reduced space
figure
hold on; grid on;
plot(Omega / (2*pi), rho, 'k', 'LineWidth', 2)
xlabel('$f$ [kHz]', 'Interpreter', 'latex')
ylabel('$\rho$', 'Interpreter', 'latex')
title('Reduced Backbone', 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
axis square tight
box on

% Plot backbone in physical space
figure
hold on; grid on;
plot(Omega / (2*pi), wBB, 'k', 'LineWidth', 2)
xlabel('$f$ [kHz]', 'Interpreter', 'latex')
ylabel('$w$ [$\mathrm{\mu}$m]', 'Interpreter', 'latex')
title('Backbone', 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
axis square tight
box on
