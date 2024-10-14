% EXAMPLE: electric field meshed with 2D element
clear
clc
close all

%% Electrostatic domain
nx = 10;
ny = 10;
Lx = 1;         % [m]
Ly = 1;         % [m]
thickness = 1;  % [m]
[nodes, elements, nset] = mesh_2Drectangle(Lx, Ly, nx, ny, 'Quad4');

% Create a dielectric material with relative permittivity of 1 [F/m]
DielMat = DielectricMaterial();
DielMat.set('RELATIVE_PERMITTIVITY', 1);

myElementConstructor = @() Quad4Element_EM(thickness, DielMat);

eMesh = Mesh(nodes);
eMesh.create_elements_table(elements, myElementConstructor);
eMesh.get_edges(1, 0);

voltage_0 = 10;
eMesh.set_essential_boundary_condition(nset{2}, 1, 0);
eMesh.set_essential_boundary_condition(nset{4}, 1, voltage_0);

% Create an assembly object for the mesh
eAssembly = Assembly(eMesh);

%% Electric potential
% Compute the electrostatic stiffness matrix
Kel = eAssembly.matrix('electrostatic_stiffness_matrix');

% Partition the stiffness matrix into submatrices based on constrained and unconstrained DOFs
Kuu = Kel(eMesh.EBC.unconstrainedDOFs, eMesh.EBC.unconstrainedDOFs);
Kuc = Kel(eMesh.EBC.unconstrainedDOFs, eMesh.EBC.constrainedDOFs(:, 1));
Kcu = Kel(eMesh.EBC.constrainedDOFs(:, 1), eMesh.EBC.unconstrainedDOFs);
Kcc = Kel(eMesh.EBC.constrainedDOFs(:, 1), eMesh.EBC.constrainedDOFs(:, 1));

% Extract the constrained voltage values
Vbc = eMesh.EBC.constrainedDOFs(:, 2);

% Solve for the electric potential at the unconstrained DOFs
potential = Kuu \ (-Kuc * Vbc);

% Unconstrain the potential vector to include constrained DOFs
potential = eAssembly.unconstrain_vector(potential);

figure(1)
PlotFieldonMesh(nodes, elements, potential, 'color', 'w');
title('Electric Potential (\phi)')
colormap hot

%% Electric field

Ex = eAssembly.vector('electric_field_x', potential);
Ey = eAssembly.vector('electric_field_y', potential);

% Ex and Ey are computed summing over the nodes, but they should be
% averaged instead. Let's count each node's appearances:
[c, ~, ic] = count_rows(sort(elements(:)));
[nodeID, ind] = unique(ic, 'stable');
count = c(ind); % how many times a node appears in the mesh?

Ex = Ex ./ count;
Ey = Ey ./ count;
E = [Ex Ey];

figure(1)
hold on
quiver(nodes(:, 1), nodes(:, 2), Ex, Ey, 0.5, 'b')
legend('\phi', 'E=-\nabla\phi', 'fontsize', 12)

figure(2)
PlotFieldonMesh(nodes, elements, sqrt(Ex.^2 + Ey.^2), 'color', 'w');
title('Electric field (E)')
colormap hot