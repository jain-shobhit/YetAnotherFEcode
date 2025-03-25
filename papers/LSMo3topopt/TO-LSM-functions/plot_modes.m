function fig = plot_modes(myAssembly, to, nev)
    % Function to plot the modes of a structure
    % Inputs:
    %   myAssembly: instance of the Assembly class
    %   to: instance of the TopologyOptimization class
    %   nev: number of modes to find
    % Outputs:
    %   fig: figure handle

    % Apply filtering and projection stages
    to.filter();
    to.projection();
    to.simp();

    % Assemble matrices
    u0 = zeros(myAssembly.Mesh.nDOFs, 1);
    [K, ~] = myAssembly.tangent_stiffness_and_force_uniform(u0, 'weights', to.d_simp);
    M = myAssembly.mass_matrix_uniform('weights', to.d_simp);

    % Constrain matrices
    Kc = myAssembly.constrain_matrix(K);
    Mc = myAssembly.constrain_matrix(M);

    % Solve eigenvalue problem
    [V, D] = eigs(Kc, Mc, nev, 'smallestabs');
    omegas = sqrt(diag(D));

    % Loop over modes
    fig = figure;
    tiledlayout(fig, 'flow', 'TileSpacing', 'compact')
    for ii = 1:nev
        iMode = ii;
        
        % Extract frequency
        omega = omegas(iMode);
        f = omega / (2 * pi);
        
        % Extract mode
        uc = V(:, iMode);
        uc = uc / sqrt(uc.' * Mc * uc);
        u = myAssembly.unconstrain_vector(uc);
        
        % Plot
        nexttile(ii);
        plot_layout_deformed(to.nel, to.d, myAssembly.Mesh.nodes, u, to.mapFea2To)
        title(num2str([iMode, f], 'Mode %d: %.2e'), 'Interpreter', 'latex')
    end
end
