function sol = mbb_gamma_topopt(myMesh, myAssembly, to, m, move, maxIter)

% Element quantities
Ke = myMesh.Elements(1).Object.tangent_stiffness_and_force(zeros(8,1));
Me = myMesh.Elements(1).Object.mass_matrix();
T2e = myMesh.Elements(1).Object.T2();
T3e = myMesh.Elements(1).Object.T3();

% Area
Ae = myMesh.Elements(1).Object.area;
Atot = Ae * myMesh.nElements;

% Null displacement vector
nDOFs = myMesh.nDOFs;
u0 = zeros(nDOFs, 1);

% History file
history = NaN(m + 1, maxIter);
densHistory = NaN(to.nElements, maxIter);

% Initialize the optimizer (MMA)
warning('off', 'all')
mma = MMA(m, move, to);

%% Main loop

% Header
fprintf("\nIteration - Frequency - Gamma - Area\n");

% Start timer
tStart = tic;

% Loop
iter = 1;
while (iter < maxIter)
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
    [idx, ~] = modal_assurance_criterion(V, modeRef);
    omega = omegas(idx);
    uc = V(:, idx);
    uc = uc / sqrt(uc.' * Mc * uc);
    u = myAssembly.unconstrain_vector(uc);
    f = omega / (2 * pi);

    % Physical density sensitivity
    sensPh = to.simp_sensitivity();

    % Compute gamma and its sensitivity
    T2 = myAssembly.tensor_uniform('T2', [nDOFs, nDOFs, nDOFs], [2, 3], 'weights', to.d_simp);
    T3 = myAssembly.tensor_uniform('T3', [nDOFs, nDOFs, nDOFs, nDOFs], [2, 3, 4], 'weights', to.d_simp);
    T2c = myAssembly.constrain_tensor(T2);
    T3c = myAssembly.constrain_tensor(T3);
    [g, dgdd] = gamma_sensitivity(myAssembly, Mc, Kc, T2c, T3c, omega, uc, Me, Ke, T2e, T3e, sensPh, sensPh);

    % Store initial values
    if iter == 1 || mod(iter, 40) == 0
        f0 = f;
        g0 = abs(g);
    end

    % Compute sensitivity
    dfdd = SensitivityLibrary.frequency(myMesh, omega, u, Ke, Me, sensPh, sensPh);
    dAdd = Ae * ones(myMesh.nElements, 1);

    % Compute filter sensitivity
    df = to.filter_sensitivity(dfdd);
    dg = to.filter_sensitivity(dgdd);
    dA = to.filter_sensitivity(dAdd);
    
    % Print current iteration
    fprintf("\n%4d %16.4e %16.4e %16.4e", iter, f, g, A / Atot);

    % Design variables (n x 1)
    xval  = to.d;

    % Objective function and sensitivity (n x 1)
    f0val = -f / f0;
    df0dx = -df(:) / f0;

    % Constraints (m x 1) and sensitivity (m x n)
    fval  = [A / Atot - 0.5; (g - (-1e-3)) / g0];
    dfdx  = [dA(:).' / Atot; dg(:).' / g0];

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

% Save optimal physical density
sol.densOpt = to.d_simp;

% Save eigenfrequency and mode shape
sol.omega = omega;
sol.phi = uc;
sol.gamma = g;

% Resize history
nPoints = find(~isnan(history(1, :)), 1, 'last');
sol.history = history(:, 1:nPoints);
sol.densHistory = densHistory(:, 1:nPoints);
end
