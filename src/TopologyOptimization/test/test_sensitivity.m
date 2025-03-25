clear; clc; close all;

%% Problem settings
% Material
E         = 70e9;   % Young's modulus [Pa]
rho       = 2700;   % density [kg/m^3]
nu        = 0.33;   % Poisson's ratio 
thickness = 0.1;    % [m] beam's out-of-plane thickness
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	% set "false" for plane_strain
myElementConstructor = @()Quad4Element(thickness, myMaterial);

% Mesh
lx = 1; ly = 0.5;
nelx = 100; nely = 50;
[nodes, elements, nset] = mesh_2Drectangle(lx, ly, nelx, nely, 'QUAD4');
myMesh = Mesh(nodes);
myMesh.create_elements_table(elements, myElementConstructor);

% Boundary conditions
myMesh.set_essential_boundary_condition(nset{1}, 1:2, 0);
myMesh.set_essential_boundary_condition(nset{3}, 1, 0);

% Assembly
myAssembly = Assembly(myMesh);

% Nodal force
F = zeros(myMesh.nDOFs,1);
nf = find_node(lx, ly, [], nodes); % node where to put the force
node_force_dofs = get_index(nf, myMesh.nDOFPerNode);
F(node_force_dofs(2)) = -1e3;
Fc = myAssembly.constrain_vector(F);

% Elements centroid
coord = zeros(myMesh.nElements, 2);
for ii = 1:myMesh.nElements
    coord(ii, :) = mean(myMesh.Elements(ii).Object.nodes);
end

% Element stiffness matrix
Ke = myMesh.Elements(1).Object.tangent_stiffness_and_force(zeros(8,1));

%% Initialize topology optimization
radius = 4;
beta = 10; eta = 0.5;
dMinSimp = 1e-6; p = 3;

% Initialize object
to = TopologyOptimization([nelx, nely], coord, radius, beta, eta, dMinSimp, p);

% Initial layout
to.initialize_density(0.5);
d0 = to.d;

% Initial layout
figure();
plot_layout(to.nel, to.d, to.mapFea2To);
title('Initial Layout', 'Interpreter', 'latex');
drawnow;

%% Analysis

% Evaluate objective function and sensitivity
[J0, dJ] = evaluate(d0, to, myAssembly, Fc, Ke, 'computeSensitivity', true);

%% Sensitivity test

% Header
fprintf("\nAnalytical, Finite Differences, Absolute Error, Relative Error\n");

% Finite differences parameters
pert = -10.^(-7:-3);
nPert = length(pert);
nEl = 50;
elements = randi(myAssembly.Mesh.nElements, 1, nEl);

% Initialize variables
errRel = zeros(nEl, nPert);

% Loop over elements
for i = 1:nEl
    % Extract element sensitivity
    el = elements(i);
    fprintf('\nElement %5d: %12.4e', el, dJ(el))
    for j = 1:nPert
        % Apply perturbation
        d = d0;
        d(el) = d(el) + pert(j);

        % Evaluate objective
        J = evaluate(d, to, myAssembly, Fc, Ke);

        % Finite differences
        sensFD = (J - J0) / pert(j);
        errAbs = abs(sensFD - dJ(el));
        errRel(i, j) = errAbs / abs(dJ(el));

        % Print
        fprintf('\n\tPerturbation %.2e: %12.4e, %12.4e, %12.4e', ...
            pert(j), sensFD, errAbs, errRel(i, j))
    end
end
fprintf('\n\nDone.\n')

% Plot
figure;
semilogy(1:nEl, errRel(:, 1), '.-', 'LineWidth', 2, 'MarkerSize', 20, ...
    'DisplayName', num2str(pert(1), '%.2e'))
hold on; grid on; box on; axis tight;
for j = 2:nPert
    semilogy(1:nEl, errRel(:, j), '.-', 'LineWidth', 2, 'MarkerSize', 20, ...
        'DisplayName', num2str(pert(j), '%.2e'))
end
xlabel('Element', 'Interpreter', 'latex')
ylabel('Relative Error', 'Interpreter', 'latex')
title('Sensitivity Test', 'Interpreter', 'latex')
legend('Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')

%% Function
function [J, dJ] = evaluate(d, to, myAssembly, Fc, Ke, varargin)
    % Parse inputs
    parser = inputParser;
    addOptional(parser, 'computeSensitivity', false);
    parse(parser, varargin{:});
    computeSensitivity = parser.Results.computeSensitivity;

    % Assign density
    to.d = d;

    % Apply filtering and projection stages
    to.filter();
    to.projection();
    to.simp();

    % Assemble matrix
    u0 = zeros(myAssembly.Mesh.nDOFs, 1);
    [K, ~] = myAssembly.tangent_stiffness_and_force_uniform(u0, 'weights', to.d_simp);
    Kc = myAssembly.constrain_matrix(K);

    % Solve stationary problem
    uc = Kc \ Fc;
    u = myAssembly.unconstrain_vector(uc);

    % Compliance
    J = dot(Fc, uc);

    % Compute sensitivity
    if computeSensitivity
        % Physical density sensitivity
        sensPh = to.simp_sensitivity();
    
        % Compute physical sensitivity
        dCdd = SensitivityLibrary.compliance(myAssembly.Mesh, u, Ke, sensPh);
    
        % Compute filter sensitivity
        dJ = to.filter_sensitivity(dCdd);
    end
end
