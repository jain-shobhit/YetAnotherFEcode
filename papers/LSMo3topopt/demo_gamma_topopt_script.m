clear; clc; close all;

addpath TO-LSM-functions
run ..\..\startup.m

%% Finite element model

% Material
E   = 148e9;    % Young's modulus [Pa]
rho = 2330e-6;  % density [ng/um^3]
nu  = 0.23;     % Poisson's ratio [-]

thickness = 24; % beam's out-of-plane thickness [um]

myMaterial = KirchoffMaterial();
set(myMaterial, 'YOUNGS_MODULUS', E, 'DENSITY', rho, 'POISSONS_RATIO', nu);
myMaterial.PLANE_STRESS = true;	% set "false" for plane_strain
myElementConstructor = @()Quad4Element(thickness, myMaterial);

% Design domain
lx = 500; ly = 100; % domain sizes [um]
nelx = 100; nely = 20;

[nodes, elements, nset] = mesh_2Drectangle(lx, ly, nelx, nely, 'QUAD4');

% Mesh
myMesh = Mesh(nodes);
myMesh.create_elements_table(elements, myElementConstructor);
nDOFs = myMesh.nDOFs;

% Boundary conditions
myMesh.set_essential_boundary_condition(nset{1}, 1:2, 0);
myMesh.set_essential_boundary_condition(nset{3}, 1, 0);

% Assembly
myAssembly = Assembly(myMesh);

% Target DOF
targetNode = find_node(lx, ly/2, [], nodes);
targetDofs = get_index(targetNode, myMesh.nDOFPerNode);
targetDof = myAssembly.free2constrained_index(targetDofs(2));

%% Topology optimization settings

% Filter parameters
radius = 4;
beta = 10; eta = 0.5;
dMin = 1e-6; p = 1;

% Elements centroid
coord = zeros(myMesh.nElements, 2);
for ii = 1:myMesh.nElements
    coord(ii, :) = mean(myMesh.Elements(ii).Object.nodes);
end

% Initialize object
to = TopologyOptimization([nelx, nely], coord, radius, beta, eta, dMin, p);

% Initial layout with uniform 0.5 density
to.initialize_density(0.5);

% Fix the proof mass with uniform 1.0 density
to.set_density_box([lx, ly], [0.2 * lx, 1e10], 1);

%% Eigenfrequency optimization

% Initialize optimizer
nConstraints = 1;
moveLimit = 0.01;
maxIter = 500;

% Optimize
solOmegaOpt = mbb_omega_topopt(myMesh, myAssembly, to, nConstraints, moveLimit, maxIter);

%% Eigenfrequency optimization results

% Optimization history
figure();
plot_history(solOmegaOpt.history);

% Optimal layout
figure();
plot_layout([nelx, nely], solOmegaOpt.densHistory(:, end), to.mapFea2To);
title('Optimal Layout', 'Interpreter', 'latex');

% Optimal density field (SIMP)
densOpt = solOmegaOpt.densOpt;

% Assemble matrices and tensors
M = myAssembly.mass_matrix_uniform('weights', densOpt);
K = myAssembly.tangent_stiffness_and_force_uniform(zeros(nDOFs, 1), 'weights', densOpt);
T2 = myAssembly.tensor_uniform('T2', [nDOFs, nDOFs, nDOFs], [2, 3], 'weights', densOpt);
T3 = myAssembly.tensor_uniform('T3', [nDOFs, nDOFs, nDOFs, nDOFs], [2, 3, 4], 'weights', densOpt);

% Apply boundary conditions
Mc = myAssembly.constrain_matrix(M);
Kc = myAssembly.constrain_matrix(K);
T2c = myAssembly.constrain_tensor(T2);
T3c = myAssembly.constrain_tensor(T3);

% Compute gamma
g = gamma_value(Mc, Kc, T2c, T3c, solOmegaOpt.omega, solOmegaOpt.phi);
fprintf('\nGamma value: %.4e\n', g)

% Compute rho from the target physical displacement
wTarget = 10;
rhoTarget = compute_rho(wTarget, Mc, Kc, T2c, T3c, solOmegaOpt.omega, solOmegaOpt.phi, targetDof);
rho = linspace(0, rhoTarget, 51);

% Compute manifold and backbone
compute_manifold(Mc, Kc, T2c, T3c, solOmegaOpt.omega, solOmegaOpt.phi, rho, targetDof, true);

%% Gamma optimization

% Initialize optimizer
nConstraints = 2;
moveLimit = 0.01;
maxIter = 500;

% Optimize
solGammaOpt = mbb_gamma_topopt(myMesh, myAssembly, to, nConstraints, moveLimit, maxIter);

%% Gamma optimization results

% Optimization history
figure();
plot_history(solGammaOpt.history);

% Optimal layout
figure();
plot_layout([nelx, nely], solGammaOpt.densHistory(:, end), to.mapFea2To);
title('Optimal Layout', 'Interpreter', 'latex');

% Optimal density field (SIMP)
densOpt = solGammaOpt.densOpt;

% Assemble matrices and tensors
M = myAssembly.mass_matrix_uniform('weights', densOpt);
K = myAssembly.tangent_stiffness_and_force_uniform(zeros(nDOFs, 1), 'weights', densOpt);
T2 = myAssembly.tensor_uniform('T2', [nDOFs, nDOFs, nDOFs], [2, 3], 'weights', densOpt);
T3 = myAssembly.tensor_uniform('T3', [nDOFs, nDOFs, nDOFs, nDOFs], [2, 3, 4], 'weights', densOpt);

% Apply boundary conditions
Mc = myAssembly.constrain_matrix(M);
Kc = myAssembly.constrain_matrix(K);
T2c = myAssembly.constrain_tensor(T2);
T3c = myAssembly.constrain_tensor(T3);

% Compute gamma
g = gamma_value(Mc, Kc, T2c, T3c, solGammaOpt.omega, solGammaOpt.phi);
fprintf('\nGamma value: %.4e\n', g)

% Compute rho from the target physical displacement
wTarget = 10;
rhoTarget = compute_rho(wTarget, Mc, Kc, T2c, T3c, solGammaOpt.omega, solGammaOpt.phi, targetDof);
rho = linspace(0, rhoTarget, 51);

% Compute manifold and backbone
compute_manifold(Mc, Kc, T2c, T3c, solGammaOpt.omega, solGammaOpt.phi, rho, targetDof, true);
