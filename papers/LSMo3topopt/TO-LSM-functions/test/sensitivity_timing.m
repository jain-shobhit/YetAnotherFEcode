% Time the evaluation for the gamma and omega sensitivities
clear; clc; close all;

% Mesh size
% nelx = 2.^(5:8);
% nelx = [40, 80, 100, 160, 200];
% nely = nelx ./ 2;
nelx = [100, 100, 100, 200];
nely = [20, 50, 100, 100];
n = nelx .* nely;

%% Sensitivity
tElapsedOmega = zeros(size(n));
tElapsedGamma = zeros(size(n));
for ii = 1:length(n)
    fprintf('\n')
    [tElapsedOmega(ii), tElapsedGamma(ii)] = evaluate_gamma_sensitivity(nelx(ii), nely(ii));
    fprintf('\tomega done in %f seconds.\n', tElapsedOmega(ii))
    fprintf('\tgamma done in %f seconds.\n', tElapsedGamma(ii))
end

%% Plot
figure
loglog(n, tElapsedOmega, '.-', 'LineWidth', 2, 'MarkerSize', 20)
hold on; grid on;
loglog(n, tElapsedGamma, '.-', 'LineWidth', 2, 'MarkerSize', 20)
xlabel('nElements', 'Interpreter', 'latex')
ylabel('Time [s]', 'Interpreter', 'latex')
legend('$\omega$', '$\gamma$', 'Location', 'Northwest', 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')

%% Evaluate gamma sensitivity
function [tElapsedOmega, tElapsedGamma] = evaluate_gamma_sensitivity(nelx, nely)
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
    % nelx = 8; nely = 4;
    [nodes, elements, nset] = mesh_2Drectangle(lx, ly, nelx, nely, 'QUAD4');
    myMesh = Mesh(nodes);
    myMesh.create_elements_table(elements, myElementConstructor);
    
    % Boundary conditions
    myMesh.set_essential_boundary_condition(nset{1}, 1:2, 0);
    myMesh.set_essential_boundary_condition(nset{3}, 1:2, 0);
    
    % Assembly
    myAssembly = Assembly(myMesh);
    nDOFs = myMesh.nDOFs;
    
    % Nodal force
    F = zeros(myMesh.nDOFs,1);
    nf = find_node(lx/2, ly/2, [], nodes); % node where to put the force
    node_force_dofs = get_index(nf, myMesh.nDOFPerNode);
    F(node_force_dofs(2)) = -1e3;
    Fc = myAssembly.constrain_vector(F);
    
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
    
    % Null displacement vector
    u0 = zeros(myMesh.nDOFs, 1);
    
    %% Initialize topology optimization
    radius = 2;
    beta = 1; eta = 0.5;
    dMin = 1e-6; p = 3;
    
    % Initialize object
    to = TopologyOptimization([nelx, nely], coord, radius, beta, eta, dMin, p);
    
    % Initial layout
    to.initialize_density(0.5);
    
    % Modify regions of the domain
    to.set_density_box([lx/2, ly/2], [0.2 * lx, 1e10], 1);
    
    %% Evaluate
    
    % Apply filtering and projection stages
    to.filter();
    to.projection();
    to.simp();

    % Physical density sensitivity
    sensPhK = to.simp_sensitivity();
    sensPhM = to.simp_sensitivity();
    
    % Assemble matrices
    [K, ~] = myAssembly.tangent_stiffness_and_force_uniform(u0, 'weights', to.d_simp);
    M = myAssembly.mass_matrix_uniform('weights', to.d_simp);    
    Kc = myAssembly.constrain_matrix(K);
    Mc = myAssembly.constrain_matrix(M);
    
    % Solve reference problem
    vc = Kc \ Fc;
    
    % Solve eigenvalue problem
    nev = 10;
    [V, D] = eigs(Mc, Kc, nev, 'LM');
    omegas = sqrt(1 ./ diag(D));
    
    % Apply MAC
    [idx, ~] = modal_assurance_criterion(V, vc);
    omega = omegas(idx);
    f = omega / (2*pi);
    uc = V(:, idx);
    uc  = uc ./ sqrt(uc.'*Mc*uc);
    u = myAssembly.unconstrain_vector(uc);

    % Compute omega
    tStart = tic;
    dfdd = SensitivityLibrary.frequency(myMesh, omega, u, Ke, Me, sensPhK, sensPhM);
    tElapsedOmega = toc(tStart);
    
    % Compute gamma
    tStart = tic;
    T2 = myAssembly.tensor_uniform('T2', [nDOFs, nDOFs, nDOFs], [2, 3], 'weights', to.d_simp);
    T3 = myAssembly.tensor_uniform('T3', [nDOFs, nDOFs, nDOFs, nDOFs], [2, 3, 4], 'weights', to.d_simp);
    T2c = myAssembly.constrain_tensor(T2);
    T3c = myAssembly.constrain_tensor(T3);
    [gamma, dgamma] = gamma_sensitivity(myAssembly, Mc, Kc, T2c, T3c, omega, uc, Me, Ke, T2e, T3e, sensPhK, sensPhM);
    tElapsedGamma = toc(tStart);
end






























