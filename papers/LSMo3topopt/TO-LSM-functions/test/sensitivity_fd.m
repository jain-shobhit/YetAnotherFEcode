% Test the sensitivity of gamma with respect to the design variables
clear; clc; close all;
addpath(genpath('..'))

%% Problem settings
% Material
E = 148e9; % Young's modulus [Pa]
rho = 2330e-6; % density [ng/um^3]
nu = 0.23; % Poisson's ratio 
thickness = 24; % [um] beam's out-of-plane thickness
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	% set "false" for plane_strain
myElementConstructor = @()Quad4Element(thickness, myMaterial);

% Mesh
lx = 500; ly = 100; % [um]
nelx = 100; nely = 20;
[nodes, elements, nset] = mesh_2Drectangle(lx, ly, nelx, nely, 'QUAD4');
myMesh = Mesh(nodes);
myMesh.create_elements_table(elements, myElementConstructor);

% Boundary conditions
myMesh.set_essential_boundary_condition(nset{1}, 1:2, 0);
myMesh.set_essential_boundary_condition(nset{3}, 1, 0);

% Assembly
myAssembly = Assembly(myMesh);
nDOFs = myMesh.nDOFs;

% Elements centroid
coord = zeros(myMesh.nElements, 2);
for ii = 1:myMesh.nElements
    coord(ii, :) = mean(myMesh.Elements(ii).Object.nodes);
end

% Element matrices
Ke = myMesh.Elements(1).Object.tangent_stiffness_and_force(zeros(8,1));
Me = myMesh.Elements(1).Object.mass_matrix();
T2e = myMesh.Elements(1).Object.T2;
T3e = myMesh.Elements(1).Object.T3;

% Area
Ae = myMesh.Elements(1).Object.area;
Atot = Ae * myMesh.nElements;

% Null displacement vector
u0 = zeros(myMesh.nDOFs, 1);

%% Initialize topology optimization
radius = 4;
beta = 10; eta = 0.5;
dMin = 1e-6; p = 3;

% Initialize object
to = TopologyOptimization([nelx, nely], coord, radius, beta, eta, dMin, p);

% Initial layout
to.initialize_density(1);

% Modify regions of the domain
to.set_initial_density_box([lx, ly], [0.2 * lx, 1e10], 1);

% Initial layout
figure();
plot_layout(to.nel, to.d, to.mapFea2To);
title('Initial Layout', 'Interpreter', 'latex');
drawnow

%% Compute sensitivity

% Initial density
d0 = to.d;

% Initial value and sensitivity
[J0, dJ, modeRef] = fun(d0, to, myAssembly, Ke, Me, T2e, T3e, 'computeSens', true);

%% Finite differences

% Perturbed elements
nIdx = 50;
elIdx = randi(myAssembly.Mesh.nElements, nIdx, 1);

% Perturbation steps
dPert = 10.^(-2:-1:-5);
nPert = length(dPert);

% Initialize
dJfd = zeros(nIdx, nPert);
errAbs = zeros(nIdx, nPert);
errRel = zeros(nIdx, nPert);

% Loop over elements
tStart = tic;
for i = 1:nIdx
    % Extract element
    thisIdx = elIdx(i);
    dJi = dJ(thisIdx);
    % Display
    fprintf('Element %d: %12.4e\n', thisIdx, dJi);

    % Loop over perturbation steps
    for j = 1:nPert
        % Negative perturbation
        d = d0;
        d(thisIdx) = d0(thisIdx) - dPert(j);

        % Compute sensitivity
        J = fun(d, to, myAssembly, Ke, Me, T2e, T3e, 'modeRef', modeRef);

        % Finite differences
        dJfd(i, j) = (J0 - J) / dPert(j);

        % Error
        errAbs(i, j) = abs(dJfd(i, j) - dJi);
        errRel(i, j) = errAbs(i, j) / abs(dJi);

        % Display
        fprintf('\tStep %.0e: %12.4e, %12.4e, %12.4e\n', dPert(j), dJfd(i, j), errAbs(i, j), errRel(i, j));
    end
end
tElapsed = toc(tStart);
fprintf('\nElapsed time is %.2f seconds.\n\n', tElapsed)

%% Plot

% Plot sensitivity
figure();
semilogy(1:length(elIdx), errRel(:, 1), '.-', 'LineWidth', 2, 'MarkerSize', 20, ...
    'DisplayName', sprintf('$\\delta\\mu = %.1e$', dPert(1)));
for j = 2:length(dPert)
    hold on; grid on; box on; axis tight;
    semilogy(1:length(elIdx), errRel(:, j), '.-', 'LineWidth', 2, 'MarkerSize', 20, ...
        'DisplayName', sprintf('$\\delta\\mu = %.1e$', dPert(j)));
end
xlabel('Element index', 'Interpreter', 'latex');
ylabel('Relative error', 'Interpreter', 'latex');
legend('Interpreter', 'latex', 'NumColumns', 2, 'Location', 'south')
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
set(gca, 'YTick', 10.^(-9:2:3))
ylim([1e-7, 1e1])

%% Save
% exportgraphics(gcf, 'sensitivity_validation.jpg', 'Resolution', 300);
% save('sensitivity_validation.mat', 'dJfd', 'errAbs', 'errRel', 'tElapsed', 'dPert', 'elIdx')

%% Function for gradient check
function [J, dJ, modeRef] = fun(x, to, myAssembly, Ke, Me, T2e, T3e, varargin)
    % Parse inputs
    p = inputParser;
    addOptional(p, 'modeRef', []);
    addOptional(p, 'computeSens', false);
    parse(p, varargin{:});
    modeRef = p.Results.modeRef;
    computeSens = p.Results.computeSens;

    % Update topology optimization
    % to.d = x;
    % to.filter();
    % to.projection();
    % to.simp();
    % mu = to.d_simp;
    % sensPhK = to.simp_sensitivity();
    % sensPhM = to.simp_sensitivity();

    % Design variables
    mu = x;
    sensPhK = ones(to.nElements, 1);
    sensPhM = ones(to.nElements, 1);

    % Number of dofs
    nDOFs = myAssembly.Mesh.nDOFs;
    u0 = zeros(nDOFs, 1);

    % Assemble matrices
    [K, ~] = myAssembly.tangent_stiffness_and_force_uniform(u0, 'weights', mu);
    M = myAssembly.mass_matrix_uniform('weights', mu);
    T2 = myAssembly.tensor_uniform('T2', [nDOFs, nDOFs, nDOFs], [2, 3], 'weights', mu);
    T3 = myAssembly.tensor_uniform('T3', [nDOFs, nDOFs, nDOFs, nDOFs], [2, 3, 4], 'weights', mu);
    
    Kc = myAssembly.constrain_matrix(K);
    Mc = myAssembly.constrain_matrix(M);
    T2c = myAssembly.constrain_tensor(T2);
    T3c = myAssembly.constrain_tensor(T3);
    
    % Solve eigenvalue problem
    nev = 10;
    [V, D] = eigs(Kc, Mc, nev, 'SM');
    omegas = sqrt(diag(D));
    
    % Check if the reference mode is provided
    if isempty(modeRef)
        modeIdx = 1;
        uc = V(:, modeIdx);
        uc  = uc ./ sqrt(uc.'*Mc*uc);
        modeRef = uc;
    else
        modeIdx = modal_assurance_criterion(V, modeRef);
    end

    % Extract target mode
    omega = omegas(modeIdx);
    uc = V(:, modeIdx);
    uc  = uc ./ sqrt(uc.'*Mc*uc);

    % Compute gamma
    if computeSens
        [J, dJdd] = gamma_sensitivity(myAssembly, Mc, Kc, T2c, T3c, omega, uc, Me, Ke, T2e, T3e, sensPhK, sensPhM);
        % dJ = to.filter_sensitivity(dJdd);
        dJ = dJdd;
    else
        J = gamma_value(Mc, Kc, T2c, T3c, omega, uc);
        dJ = [];
    end
end
