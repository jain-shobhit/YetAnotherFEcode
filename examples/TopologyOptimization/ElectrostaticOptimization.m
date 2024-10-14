clear; clc; close all;

%% Problem settings
% Material
thickness = 1;
DielMat = DielectricMaterial();
DielMat.set('RELATIVE_PERMITTIVITY', 1);
myElementConstructor = @() Quad4Element_EM(thickness, DielMat);

% Mesh
lx = 1; ly = 1;
nelx = 100; nely = 100;
[nodes, elements, nset] = mesh_2Drectangle(lx, ly, nelx, nely, 'QUAD4');
myMesh = Mesh(nodes);
myMesh.create_elements_table(elements, myElementConstructor);

% Boundary conditions
nodesTemp = 1:myMesh.nNodes;
boundaryNodes = nodesTemp(all(abs(myMesh.nodes - [0.5 * lx, 0]) < [0.1 * lx, 1e-10], 2));
boundaryDofs = get_index(boundaryNodes, myMesh.nDOFPerNode);
myMesh.set_essential_boundary_condition(boundaryDofs, 1, 0);

% Assembly
myAssembly = Assembly(myMesh);

% Load
F = 1e-6 * ones(myMesh.nDOFs, 1);
Fc = myAssembly.constrain_vector(F);

% Elements centroid
coord = zeros(myMesh.nElements, 2);
for ii = 1:myMesh.nElements
    coord(ii, :) = mean(myMesh.Elements(ii).Object.nodes);
end

% Element stiffness matrix
Ke = myMesh.Elements(1).Object.electrostatic_stiffness_matrix();

% Area
Ae = myMesh.Elements(1).Object.area;
% Ae = (lx * ly) / (nelx * nely);
Atot = Ae * myMesh.nElements;

%% Initialize topology optimization
radius = 2;
beta = 1; eta = 0.5;
dMinSimp = 1e-7; p = 3;

% Initialize object
to = TopologyOptimization([nelx, nely], coord, radius, beta, eta, dMinSimp, p);

% Initial layout
to = to.initialize_density(0.5);

% Initial layout
figure();
plot_layout(to.nel, to.d, to.mapFea2To);
title('Initial Layout', 'Interpreter', 'latex');
drawnow;

%% Initialize optimizer
m = 1;
move = 0.01;
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
tic;

% Loop
iter = 1;
while (iter < maxIter)
    % Update beta
    if mod(iter, 50) == 0
        to.beta = to.beta * 1.5;
        fprintf("\nUpdate beta: %.4e", to.beta);
    end

    % Apply filtering and projection stages
    to = to.filter();
    to = to.projection();
    to = to.simp();

    % Current area
    A = Ae * sum(to.d_proj);

    % Assemble matrix
    K = myAssembly.matrix_uniform('electrostatic_stiffness_matrix', 'weights', to.d_simp);
    Kc = myAssembly.constrain_matrix(K);

    % Solve stationary problem
    uc = Kc \ Fc;
    u = myAssembly.unconstrain_vector(uc);
    C = dot(Fc, uc);

    % Store initial value
    if iter == 1
        C0 = C;
    end

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
    f0val = C / C0;
    df0dx = dC(:) / C0;

    % Constraints (m x 1) and sensitivity (m x n)
    fval  = [A / Atot - 0.5];
    dfdx  = [dA(:).' / Atot];

    % Save current iteration
    history(:, iter) = [f0val, fval];
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
    [xmma, mma] = mma.optimize(iter, xval, f0val, df0dx, fval, dfdx);
    to.d = xmma;

    % Update counter
    iter = iter + 1;
end

% Stop timer and display elapsed time
fprintf('\n\nEnd of the optimization.\n');
toc;

%% Optimal results

% History
figure();
plot_history(history);

% Optimal layout
figure();
plot_layout(to.nel, to.d_proj, to.mapFea2To);
title('Optimal Layout', 'Interpreter', 'latex');
