function myAssembly = define_beam_model(p)
% Define the model for the von Karman beam using Yet Another FE Code (YAFEC).
% The model is parametrized by the vector p.
% p = [v1, v2, h, L]
% v1: amplitude of the first mode (half-sine wave)
% v2: amplitude of the second mode (sine wave)
% h: thickness of the beam
% L: length of the beam
% Inputs:
%   p: the design variables.
% Outputs:
%   myAssembly: the assembly object.

% Geometry
L = p(4); % length
h = p(3); % in-plane thickness
b = 1e-3; % out-of-plane thickness

% mesh
nElements = 10;
dx = L/nElements;

% Nodal coordinates (nominal)
x0 = (0:dx:L).';
y0 = zeros(size(x0));

% Parametrized shape
v1 = p(1) * sin(pi / L * x0);
v2 = p(2) * sin(2 * pi / L * x0);

% Nodal coordinates (deformed)
x = x0;
y = y0 + v1 + v2;

% Material properties
E = 90e9; % Young's modulus
rho = 7850; % density
nu = 0.3; % Poisson's ratio

% Define the material
myMaterial  = KirchoffMaterial();
set(myMaterial, 'YOUNGS_MODULUS', E, 'DENSITY', rho, 'POISSONS_RATIO', nu);

% Element (same element all across the domain)
myElementConstructor = @()BeamElement(b, h, myMaterial);

nNodes = size(x,1);
nodes = [x, y];
elements = [1:nNodes-1;2:nNodes].';
beamMesh = Mesh(nodes);
beamMesh.create_elements_table(elements, myElementConstructor);

% Boundary conditions (clamped-clamped)
dirichletNodes = [1 nNodes];
constrainedDof = [1 2 3]; % u,v,theta
value = 0;
beamMesh.set_essential_boundary_condition(dirichletNodes, constrainedDof, value)

% Assembly
myAssembly = Assembly(beamMesh, false);
myAssembly.DATA.L = L;
myAssembly.DATA.b = b;
myAssembly.DATA.h = h;
