% stiffness_matrix_derivative
%
% Synthax:
% dKdq = stiffness_matrix_derivative(myAssembly, elements, V)
%
% Description: this is a wrapper function for YetAnotherFEcode that 
% computes the derivative of the tangent stiffness matrix with respect to
% the amplitude of the vibration mode V.
% INPUTS
%   - myAssembly: Assembly from YetAnotherFEcode. It MUST contain the field
%     ".DATA.K" storing the unconstrained linear stiffness matrix K0
%   - elements: table of the elements
%   - V: vector with the displacement field corresponding to a selected
%     vibration mode (V must be unconstrained)
% OUTPUTS
%   - dKdq: tangent stiffness matrix derivative wrt q, being q the
%     modal coordinate of the vibration mode V
%
% Additional notes:
%   - ALL the elements are assumed to have the same properties in terms
%     of MATERIAL and QUADRATURE rules.
%   - this function uses the stiffness_matrix_derivative function 
%     implemented in the Julia module "DpROM.jl". 
%   - as such, this function supports ONLY models meshed with the elements
%     supported by both YetAnotherFEcode AND the DpROM.jl
%   - List of currently supported elements: 
%     Q8, TET10, HEX20, WED15               (in YetAnotherFEcode)
%     Q8, TET10, HEX20, WED15, Q4, HEX8     (in DpROM.jl)
%
% Reference: J. Marconi, P. Tiso, D.E. Quadrelli & F. Braghin, "A higher 
% order parametric nonlinear reduced order model for imperfect structures 
% using Neumann expansion", Nonlinear Dynamics, 2021.
%
% Created: 14 May 2021
% Author: Jacopo Marconi, Politecnico di Milano

function dKdq = stiffness_matrix_derivative(myAssembly, elements, V)

t0=tic;

% data from myAssembly
nodes    = myAssembly.Mesh.nodes;                   % nodes table
nel      = myAssembly.Mesh.nElements;           	% number of elements
nnodes   = myAssembly.Mesh.nNodes;               	% number of nodes
freedofs = myAssembly.Mesh.EBC.unconstrainedDOFs;   % free DOFs

% Element data (assuming ALL the elements have the same properties in terms
% of material and quadrature rule)
Element = myAssembly.Mesh.Elements(1).Object;   % first element
XGauss = Element.quadrature.X;                  % quadrature points
WGauss = Element.quadrature.W;                  % quadrature weights
C = Element.initialization.C;                   % constitutive matrix

DIM = size(nodes,2);                            % 2D/3D problem
ndofs = nnodes * DIM;                           % total number of DOFs
nfree = length(freedofs);                       % number of free DOFs

% build connectivity table
conn = zeros(nel, length(Element.iDOFs));
for ee = 1 : nel
    conn(ee, :) = myAssembly.Mesh.Elements(ee).Object.iDOFs;
end

elements = int64(elements); % convert into integers
conn = int64(conn);         % convert into integers

% INPUT SIZE CHECK ________________________________________________________
% for tensor assembly, always use the full DOFs displacement vectors
if size(V,1)==nfree
    U = zeros(ndofs,size(V,2));
    U(freedofs,:) = V;
    V = U;
elseif size(V,1)~=ndofs
    error('Wrong dimensions for V')
end

% JULIA ___________________________________________________________________
% add current path in Julia
filename = 'wrapper\stiffness_matrix_derivative';   % this function name (up to folder "DpROM")
filepath = mfilename('fullpath');                   % full path to this file
folderpath = filepath(1:end-length(filename)-1);    % path to folder where "DpROM.jl" is
a = jl.eval('LOAD_PATH');                           % lists of paths already in Julia
% check if path is already present
path_already_present = 0;
for ii = 1 : length(a)
    if strcmpi(a{ii}, folderpath)==1
        path_already_present = 1;
        break
    end
end
% if not, add it
if path_already_present==0
    fp = strrep(folderpath,'\','\\'); % change separator
    jl.eval(['push!(LOAD_PATH, "' fp '")']);
    fprintf(' Path added\n\n')
end
% load packages
jleval using TensorOperations
jleval using DpROM

fprintf(' dKdq, assembling %d elements ...', nel)

% call the function in Julia to compute all the tensors
dKdq = jl.call('stiffness_matrix_derivative', elements, ...
    nodes, conn, C, V, XGauss, WGauss);

dKdq = sparse(dKdq);

fprintf(' %.2f s\n',toc(t0))



