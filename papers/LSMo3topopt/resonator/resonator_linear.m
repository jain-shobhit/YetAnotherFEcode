clear; clc; close all;
addpath(genpath(fullfile('..', 'TO-LSM-functions')))

% MEMS units:
% time [ms]
% length [um]
% mass [ng]
% frequency [kHz]
% pressure [Pa]
% density [ng/um^3]
% force [pN]

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
lx = 1000; ly = 500; % [um]
nelx = 200; nely = 100;
[nodes, elements, nset] = mesh_2Drectangle(lx, ly, nelx, nely, 'QUAD4');
myMesh = Mesh(nodes);
myMesh.create_elements_table(elements, myElementConstructor);

% Boundary conditions
myMesh.set_essential_boundary_condition(nset{1}, 1:2, 0);
myMesh.set_essential_boundary_condition(nset{3}, 1:2, 0);

% Assembly
myAssembly = Assembly(myMesh);
nDOFs = myMesh.nDOFs;

% Load
F = zeros(myMesh.nDOFs, 1);
loadNode = find_node(0.5 * lx, 0.5 * ly, [], nodes); % bottom mass
loadDofs = get_index(loadNode, myMesh.nDOFPerNode);
F0 = 1e12;
F(loadDofs(2)) = F0;
Fc = myAssembly.constrain_vector(F);

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
dMinRamp = 1e-6; q = 4;

% Initialize object
to = TopologyOptimization([nelx, nely], coord, radius, beta, eta, dMinSimp, p, dMinRamp, q);

% Initial layout
to.initialize_density(0.5);

% Modify regions of the domain
to.set_density_box([lx/2, ly/2], [0.2 * lx, 1e10], 1);

% Initialize symmetry object
symmetry_map = SymmetryMap(to.coord, to.nel, 'xy');

% Initial layout
figure();
plot_layout(to.nel, to.d, to.mapFea2To);
title('Initial Layout', 'Interpreter', 'latex');
drawnow

%% Initialize optimizer
m = 3;
move = 0.01;
% mma = MMA(m, move, to);
mma = MMA(m, move, to, symmetry_map);

% Iterations
maxIter = 300;

% History file
history = NaN(m + 1, maxIter);
densHistory = NaN(to.nElements, maxIter);

% Initialize figure
figure(); drawnow;

%% Initial modal analysis
% fig = plot_modes(myAssembly, to, 20);
% close(fig);

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
    % if mod(iter, 40) == 0
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

    % Physical density sensitivity
    sensPhK = to.simp_sensitivity();
    sensPhM = to.simp_sensitivity();

    % Current area
    A = Ae * sum(to.d_proj);

    % Assemble matrices
    [K, ~] = myAssembly.tangent_stiffness_and_force_uniform(u0, 'weights', to.d_simp);
    M = myAssembly.mass_matrix_uniform('weights', to.d_simp);

    Kc = myAssembly.constrain_matrix(K);
    Mc = myAssembly.constrain_matrix(M);

    % Solve stationary problem
    uc = Kc \ Fc;
    u = myAssembly.unconstrain_vector(uc);
    C = dot(uc, Fc) / F0;

    % Solve eigenvalue problem
    nev = 20;
    [V, D] = eigs(Kc, Mc, nev, 'smallestabs');
    omegas = sqrt(diag(D));

    % Identify mode
    if iter == 1
        modeX = V(:, 2);
        modeY = V(:, 1);
    end
    
    % Apply MAC for mode X
    [idxX, macX] = modal_assurance_criterion(V, modeX);
    omegaX = omegas(idxX);
    uXc = V(:, idxX);
    uXc = uXc / sqrt(uXc.' * Mc * uXc);
    uX = myAssembly.unconstrain_vector(uXc);
    fX = omegaX / (2 * pi);
    
    % Apply MAC for mode Y
    [idxY, macY] = modal_assurance_criterion(V, modeY);
    omegaY = omegas(idxY);
    uYc = V(:, idxY);
    uYc = uYc / sqrt(uYc.' * Mc * uYc);
    uY = myAssembly.unconstrain_vector(uYc);
    fY = omegaY / (2 * pi);

    % Store initial values
    if iter == 1 || mod(iter, 40) == 0
        C0 = C;
        fX0 = fX;
        fY0 = fY;
    end

    % Compute sensitivity
    % dCdd = SensitivityLibrary.compliance(myMesh, u, Ke, sensPhK) / F0;
    dfXdd = SensitivityLibrary.frequency(myMesh, omegaX, uX, Ke, Me, sensPhK, sensPhM);
    dfYdd = SensitivityLibrary.frequency(myMesh, omegaY, uY, Ke, Me, sensPhK, sensPhM);
    dAdd = Ae * ones(myMesh.nElements, 1);

    % Compute filter sensitivity
    % dC = to.filter_sensitivity(dCdd);
    dfX = to.filter_sensitivity(dfXdd);
    dfY = to.filter_sensitivity(dfYdd);
    dA = to.filter_sensitivity(dAdd);
    
    % Print current iteration
    fprintf("\n%4d %16.4e %16.4e %16.4e %16.4e", iter, C, fX, fY, A / Atot);

    % Plot current layout
    plot_layout(to.nel, to.d_proj, to.mapFea2To); drawnow;

    % Design variables (n x 1)
    xval  = to.d;

    % Objective function and sensitivity (n x 1)
    f0val = -fX / fX0;
    df0dx = -dfX(:) / fX0;
    % f0val = C / C0;
    % df0dx = dC(:) / C0;

    % Constraints (m x 1) and sensitivity (m x n)
    % fval  = [(C - 5) / C0; (fY / 1000 - 1).^2 - 1e-4; A / Atot - 0.6];
    % dfdx  = [dC(:).' / C0; 2 * (fY / 1000 - 1) / 1000 * dfY(:).'; dA(:).' / Atot];
    fval  = [(fY - 1010) / fY0; ...
             (990 - fY) / fY0; ...
             A / Atot - 0.6];
    dfdx  = [dfY(:).' / fY0; ...
             -dfY(:).' / fY0; ...
             dA(:).' / Atot];

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

%% Plot modes

% Create figure
figure
tiledlayout('flow', 'TileSpacing', 'compact')

% Loop over modes
modeIdx = [idxX, idxY];
for ii = 1:length(modeIdx)
    iMode = modeIdx(ii);

    % Extract frequency
    omega = omegas(iMode);
    f = omega / (2 * pi);

    % Extract mode
    uc = V(:, iMode);
    uc = uc / sqrt(uc.' * Mc * uc);
    u = myAssembly.unconstrain_vector(uc);

    % Plot
    nexttile
    plot_layout_deformed(to.nel, to.d, myMesh.nodes, u, to.mapFea2To)
    title(num2str([iMode, f], 'Mode %d: %.2e kHz'), 'Interpreter', 'latex')
end

%% Nonlinear analysis

% Evaluate gamma
fprintf('\nComputing gamma ... ')
T2 = myAssembly.tensor_uniform('T2', [nDOFs, nDOFs, nDOFs], [2, 3], 'weights', to.d_simp);
T3 = myAssembly.tensor_uniform('T3', [nDOFs, nDOFs, nDOFs, nDOFs], [2, 3, 4], 'weights', to.d_simp);
T2c = myAssembly.constrain_tensor(T2);
T3c = myAssembly.constrain_tensor(T3);
gX = gamma_value(Mc, Kc, T2c, T3c, omegaX, uXc);
gY = gamma_value(Mc, Kc, T2c, T3c, omegaY, uYc);
fprintf('done.\nGamma X: %.4e\nGamma Y: %.4e\n\n', gX, gY)
