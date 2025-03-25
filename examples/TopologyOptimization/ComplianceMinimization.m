clear; clc; close all;

%% Problem settings
% Material
E       = 70e9;     % Young's modulus [Pa]
rho     = 2700;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio 
thickness = .1;     % [m] beam's out-of-plane thickness
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

% Area
Ae = myMesh.Elements(1).Object.area;
Atot = Ae * myMesh.nElements;

% Null displacement vector
u0 = zeros(myMesh.nDOFs, 1);

%% Initialize topology optimization
radius = 2;
beta = 10; eta = 0.5;
dMinSimp = 1e-7; p = 3;

% Initialize object
to = TopologyOptimization([nelx, nely], coord, radius, beta, eta, dMinSimp, p);

% Initial layout
to.initialize_density(0.5);

% Initial layout
figure();
plot_layout(to.nel, to.d, to.mapFea2To);
title('Initial Layout', 'Interpreter', 'latex');
drawnow;

%% Initialize optimizer
m = 1;
move = 0.2;
mma = MMA(m, move, to);

% Iterations
maxIter = 200;

% History file
history = NaN(m + 1, maxIter);
densHistory = NaN(to.nElements, maxIter);

% Initialize figure
figure(); drawnow;

%% Main loop

% Header
fprintf("\nIteration - Objective - Constraints\n");

% Start timer
tStart = tic;

% Loop
iter = 1;
while (iter < maxIter)
    % Apply filtering and projection stages
    to.filter();
    to.projection();
    to.simp();

    % Current area
    A = Ae * sum(to.d_proj);

    % Assemble matrix
    [K, ~] = myAssembly.tangent_stiffness_and_force_uniform(u0, 'weights', to.d_simp);
    Kc = myAssembly.constrain_matrix(K);

    % Solve stationary problem
    uc = Kc \ Fc;
    u = myAssembly.unconstrain_vector(uc);

    % Compliance
    C = dot(Fc, uc);

    % Physical density sensitivity
    sensPh = to.simp_sensitivity();

    % Compute physical sensitivity
    dCdd = SensitivityLibrary.compliance(myMesh, u, Ke, sensPh);
    dAdd = Ae * ones(myMesh.nElements, 1);

    % Compute filter sensitivity
    dC = to.filter_sensitivity(dCdd);
    dA = to.filter_sensitivity(dAdd);
    
    % Print current iteration
    fprintf("\n%4d %16.4e %16.4e", iter, C, A / Atot);

    % Plot current layout
    plot_layout(to.nel, to.d_proj, to.mapFea2To); drawnow;

    % Design variables (n x 1)
    xval  = to.d;

    % Objective function and sensitivity (n x 1)
    f0val = C;
    df0dx = dC(:);

    % Constraints (m x 1) and sensitivity (m x n)
    fval  = [A - 0.5 * Atot];
    dfdx  = [dA(:).'];

    % Save current iteration
    history(:, iter) = [f0val; fval];
    densHistory(:, iter) = to.d_proj;
    
    % Convergence criterion
    if iter > 5 % check convergence after 5 iterations
        fval_tol = 1e-3;
        if all(fval < fval_tol) % all the constraints are satisfied
            err = abs(1 - history(1, iter-3:iter-1) / history(1, iter));
            err_tol = 1e-3;
            if all(err < err_tol) % check relative variation
                break
            end
        end
    end

    % MMA step
    to.d = mma.optimize(iter, xval, f0val, df0dx, fval, dfdx);

    % Update counter
    iter = iter + 1;
end

% Stop timer and display elapsed time
tElapsed = toc(tStart);
fprintf('\n\nEnd of the optimization.\n');
fprintf('Elapsed time is %f seconds.\n', tElapsed)

%% Optimal results

% History
figure();
plot_history(history);

% Optimal layout
figure();
plot_layout(to.nel, to.d_proj, to.mapFea2To);
title('Optimal Layout', 'Interpreter', 'latex');

% Create gif of the density evolution
create_gif(to.nel, densHistory, 'mapFea2To', to.mapFea2To, 'fileName', 'ComplianceMaximization');
