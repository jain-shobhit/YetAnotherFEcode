clear; clc; close all;
addpath(genpath(fullfile('..', 'TO-LSM-functions')))

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

% Element quantities
Ke = myMesh.Elements(1).Object.tangent_stiffness_and_force(zeros(8,1));
Me = myMesh.Elements(1).Object.mass_matrix();
T2e = myMesh.Elements(1).Object.T2();
T3e = myMesh.Elements(1).Object.T3();

% Area
Ae = myMesh.Elements(1).Object.area;
Atot = Ae * myMesh.nElements;

% Null displacement vector
u0 = zeros(myMesh.nDOFs, 1);

%% Initialize topology optimization
radius = 4;
beta = 10; eta = 0.5;
dMinSimp = 1e-6; p = 1;

% Initialize object
to = TopologyOptimization([nelx, nely], coord, radius, beta, eta, dMinSimp, p);

% Initial layout
to.initialize_density(0.5);

% Modify regions of the domain
to.set_density_box([lx, ly], [0.2 * lx, 1e10], 1);

% Initial layout
figure();
plot_layout(to.nel, to.d, to.mapFea2To);
title('Initial Layout', 'Interpreter', 'latex');
drawnow

%% Initialize optimizer
m = 2;
move = 0.01;
mma = MMA(m, move, to);

% Iterations
maxIter = 500;

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
    tIter = tic;
    % Update beta
    % if mod(iter, 50) == 0
    %     to.beta = to.beta * 1.5;
    %     fprintf("\nUpdate beta: %.4e", to.beta);
    % end
    % Update p
    if mod(iter, 40) == 0 && to.p < 8
        to.p = to.p + 1;
        fprintf("\nUpdate p: %.4e", to.p);
    end

    % Apply filtering and projection stages
    to.filter();
    to.projection();
    to.simp();

    % Current area
    A = Ae * sum(to.d_proj);

    % Assemble matrices
    [K, ~] = myAssembly.tangent_stiffness_and_force_uniform(u0, 'weights', to.d_simp);
    M = myAssembly.mass_matrix_uniform('weights', to.d_simp);

    Kc = myAssembly.constrain_matrix(K);
    Mc = myAssembly.constrain_matrix(M);

    % Solve eigenvalue problem
    nev = 10;
    [V, D] = eigs(Kc, Mc, nev, 'smallestabs');
    omegas = sqrt(diag(D));

    % Find reference mode shape
    if iter == 1
        modeRef = V(:, 1);
        modeRef = modeRef / sqrt(modeRef.' * Mc * modeRef);
    end
    
    % Apply MAC
    [idx, mac] = modal_assurance_criterion(V, modeRef);
    omega = omegas(idx);
    uc = V(:, idx);
    uc = uc / sqrt(uc.' * Mc * uc);
    u = myAssembly.unconstrain_vector(uc);
    f = omega / (2 * pi);

    % Physical density sensitivity
    sensPhK = to.simp_sensitivity();
    sensPhM = to.simp_sensitivity();

    % Compute gamma and its sensitivity
    T2 = myAssembly.tensor_uniform('T2', [nDOFs, nDOFs, nDOFs], [2, 3], 'weights', to.d_simp);
    T3 = myAssembly.tensor_uniform('T3', [nDOFs, nDOFs, nDOFs, nDOFs], [2, 3, 4], 'weights', to.d_simp);
    T2c = myAssembly.constrain_tensor(T2);
    T3c = myAssembly.constrain_tensor(T3);
    [g, dgdd] = gamma_sensitivity(myAssembly, Mc, Kc, T2c, T3c, omega, uc, Me, Ke, T2e, T3e, sensPhK, sensPhM);

    % Store initial values
    if iter == 1 || mod(iter, 40) == 0
        f0 = f;
        g0 = abs(g);
    end

    % Compute sensitivity
    dfdd = SensitivityLibrary.frequency(myMesh, omega, u, Ke, Me, sensPhK, sensPhM);
    dAdd = Ae * ones(myMesh.nElements, 1);

    % Compute filter sensitivity
    df = to.filter_sensitivity(dfdd);
    dg = to.filter_sensitivity(dgdd);
    dA = to.filter_sensitivity(dAdd);
    
    % Print current iteration
    fprintf("\n%4d %16.4e %16.4e %16.4e", iter, f, g, A / Atot);

    % Plot current layout
    plot_layout(to.nel, to.d_proj, to.mapFea2To); drawnow;

    % Design variables (n x 1)
    xval  = to.d;

    % Objective function and sensitivity (n x 1)
    f0val = -f / f0;
    df0dx = -df(:) / f0;

    % Hardening constraints (m x 1) and sensitivity (m x n)
    fval  = [(1e-3 - g) / g0; A / Atot - 0.5];
    dfdx  = [-dg(:).' / g0; dA(:).' / Atot];

    % Softening constraints (m x 1) and sensitivity (m x n)
    % fval  = [(g - (-1e-3)) / g0; A / Atot - 0.5];
    % dfdx  = [dg(:).' / g0; dA(:).' / Atot];

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
    fprintf('\b\t')
    toc(tIter)
    fprintf('\b')
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
