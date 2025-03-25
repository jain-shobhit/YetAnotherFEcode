% Test the sensitivity of gamma with respect to the design variables
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
nelx = 8; nely = 4;
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

% Area
Ae = myMesh.Elements(1).Object.area;
Atot = Ae * myMesh.nElements;

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
to.set_density_box([lx/2, ly/2], [0.2 * lx, 1e10], 0.9);

% Initialize symmetry object
symmetry_map = SymmetryMap(to.coord, to.nel);

% Initial layout
figure();
plot_layout(to.nel, to.d, to.mapFea2To);
title('Initial Layout');
drawnow

%% Main loop

% Apply filtering and projection stages
to.filter();
to.projection();
to.simp();

% Current area
A = Ae * sum(to.d_proj);

% Check gradients
opts = optimoptions('fmincon', 'FiniteDifferenceStepSize', 1e-6);
valid = checkGradients(@(x) fun(x, myAssembly, u0, Fc), to.d_simp, ...
    'Display', 'on', 'Tolerance', 1e-4);

%% Function for gradient check
function [fval, fgrad] = fun(x, myAssembly, u0, Fc)
    % Number of dofs
    nDOFs = myAssembly.Mesh.nDOFs;

    % Element matrices
    Ke = myAssembly.Mesh.Elements(1).Object.tangent_stiffness_and_force(zeros(8,1));
    Me = myAssembly.Mesh.Elements(1).Object.mass_matrix();
    T2e = myAssembly.Mesh.Elements(1).Object.T2;
    T3e = myAssembly.Mesh.Elements(1).Object.T3;

    % Assemble matrices
    [K, ~] = myAssembly.tangent_stiffness_and_force_uniform(u0, 'weights', x);
    M = myAssembly.mass_matrix_uniform('weights', x);
    T2 = myAssembly.tensor_uniform('T2', [nDOFs, nDOFs, nDOFs], [2, 3], 'weights', x);
    T3 = myAssembly.tensor_uniform('T3', [nDOFs, nDOFs, nDOFs, nDOFs], [2, 3, 4], 'weights', x);
    
    Kc = myAssembly.constrain_matrix(K);
    Mc = myAssembly.constrain_matrix(M);
    T2c = myAssembly.constrain_tensor(T2);
    T3c = myAssembly.constrain_tensor(T3);
    
    % Solve reference problem
    vc = Kc \ Fc;
    
    % Solve eigenvalue problem
    nev = 10;
    [V, D] = eigs(Mc, Kc, nev, 'LM');
    omegas = sqrt(1 ./ diag(D));
    
    % Apply MAC
    [idx, ~] = modal_assurance_criterion(V, vc);
    omega = omegas(idx);
    % f = omega / (2*pi);
    uc = V(:, idx);
    uc  = uc ./ sqrt(uc.'*Mc*uc);
    % u = myAssembly.unconstrain_vector(uc);
    
    % Compute gamma
    [gamma, dgamma] = gamma_sensitivity(myAssembly, Mc, Kc, T2c, T3c, omega, uc, Me, Ke, T2e, T3e);

    % Return
    fval = gamma;
    fgrad = dgamma;
end
