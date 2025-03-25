clear; clc; close all;

%% Problem settings
% Material
myMaterial = DielectricMaterial();
myMaterial.set('RELATIVE_PERMITTIVITY', 1e9);
myElementConstructor = @()Hex8Element_EM(myMaterial);

% Mesh
lx = 1.0; ly = 1.0; lz = 1.0;
nx = 30; ny = 30; nz = 30;
[nodes, elements, nset] = mesh_3Dparallelepiped('HEX8', lx, ly, lz, nx, ny, nz);
myMesh = Mesh(nodes);
myMesh.create_elements_table(elements, myElementConstructor);

% Boundary conditions
nodesTemp = 1:myMesh.nNodes;
boundaryNodes = nodesTemp(all(abs(myMesh.nodes - [lx/2, ly/2, 0]) < [lx/10, ly/10, 1e-10], 2));
boundaryDofs = get_index(boundaryNodes, myMesh.nDOFPerNode);
myMesh.set_essential_boundary_condition(boundaryDofs, 1, 0);

% Assembly
myAssembly = Assembly(myMesh);

% Nodal force
F = 1e-5 * ones(myMesh.nDOFs, 1);
Fc = myAssembly.constrain_vector(F);

% Elements centroid
coord = zeros(myMesh.nElements, myMesh.nDim);
for ii = 1:myMesh.nElements
    coord(ii, :) = mean(myMesh.Elements(ii).Object.nodes);
end

% Element stiffness matrix
Ke = myMesh.Elements(1).Object.electrostatic_stiffness_matrix();

% Area
Ve = myMesh.Elements(1).Object.vol;
Vtot = Ve * myMesh.nElements;

%% Initialize topology optimization
radius = 2;
beta = 10; eta = 0.5;
dMinSimp = 1e-6; p = 3;

% Initialize object
to = TopologyOptimization([nx, ny, nz], coord, radius, beta, eta, dMinSimp, p);

% Initial layout
to.initialize_density(0.5);

% Set density
to.set_density_box([lx/2, ly/2, 0], [lx/10, ly/10, lz/10], 1);
to.set_density_box([lx/2, ly/2, lz/2], lx/10, 0);

% Initial layout
figure();
plot_layout_3d([lx, ly, lz], to.coord, to.d)
drawnow;

%% Initialize optimizer
m = 1;
move = 0.05;
mma = MMA(m, move, to);

% Iterations
maxIter = 100;

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
    % Start iteration timer
    tIter = tic;

    % Apply filtering and projection stages
    to.filter();
    to.projection();
    to.simp();

    % Current area
    V = Ve * sum(to.d_proj);

    % Assemble matrix
    K = myAssembly.matrix_uniform('electrostatic_stiffness_matrix', 'weights', to.d_simp);
    Kc = myAssembly.constrain_matrix(K);

    % Solve stationary problem
    uc = Kc \ Fc;
    u = myAssembly.unconstrain_vector(uc);

    % Compliance
    C = dot(Fc, uc);
    if iter == 1
        C0 = C;
    end

    % Physical density sensitivity
    sensPh = to.simp_sensitivity();

    % Compute physical sensitivity
    dCdd = SensitivityLibrary.compliance(myMesh, u, Ke, sensPh);
    dVdd = Ve * ones(myMesh.nElements, 1);

    % Compute filter sensitivity
    dC = to.filter_sensitivity(dCdd);
    dV = to.filter_sensitivity(dVdd);
    
    % Print current iteration
    fprintf("\n%4d %16.4e %16.4e", iter, C / C0, V / Vtot);

    % Plot current layout
    plot_layout_3d([lx, ly, lz], to.coord, to.d_simp)
    drawnow;

    % Design variables (n x 1)
    xval  = to.d;

    % Objective function and sensitivity (n x 1)
    f0val = C / C0;
    df0dx = dC(:) / C0;

    % Constraints (m x 1) and sensitivity (m x n)
    fval  = V / Vtot - 0.5;
    dfdx  = dV(:).' / Vtot;

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

    % Stop iteration timer and display elapsed time
    tElapsedIter = toc(tIter);
    fprintf('\tElapsed time is %f seconds.', tElapsedIter)
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
plot_layout_3d([lx, ly, lz], to.coord, to.d_simp)
title('Optimal Layout', 'Interpreter', 'latex');
