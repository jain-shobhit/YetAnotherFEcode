% reduced_tensors_DpROM
%
% Synthax:
% tensors = reduced_tensors_DpROM(myAssembly, elements, V, U, ...
%                                 FORMULATION, volume, USEJULIA)
%
% Description: this is a wrapper function for YetAnotherFEcode that
% computes the parametric reduced order tensors for shape defect described
% in the paper.
% INPUTS
%   - myAssembly: Assembly from YetAnotherFEcode.
%   - elements: table of the elements
%   - V: Reduced Order Basis (unconstrained)
%   - U: tall matrix collecting by columns the (unconstrained) defect shapes
%   - formulation (optional): choose the Neumann formulation between "N1",
%     "N1T" and "N0".
%     The default value is N1.
%   - volume (optional): set to 1 to compute the additional tensors required
%     to approximate the volume integral over the defected volume, set to 0
%     if the integral can be computed over the nominal volume (e.g. 
%     isochoric defect).
%     The default value is 1.
%   - USEJULIA (optional, default = 0): use DpROM.jl instead of matlab
% OUTPUT
%   tensors: a struct variable with the following fields*:
%       .Q2n            n*n         reduced stiffness tensor
%    	.Q3n            n*n*n       reduced stiffness tensor
%   	.Q4n            n*n*n*n     reduced stiffness tensor
%    	.Q3d            n*n*d       reduced stiffness tensor
%       .Q4d            n*n*n*d     reduced stiffness tensor
%       .Q5d            n*n*n*n*d   reduced stiffness tensor
%      	.Q4dd           n*n*d*d     reduced stiffness tensor
%      	.Q5dd           n*n*n*d*d   reduced stiffness tensor
%     	.Q6dd           n*n*n*n*d*d reduced stiffness tensor
%       .M              n*n         reduced mass tensors**
%      	.time           computational time
%      	.formulation    adopted formulation
%      	.volume         1/0
%   *being n=size(V,2) and d=size(U,2)
%   ** parametrized mass matrix, computed using an approximated integration
%      over the defected volume. The defected mass matrix can be obtained 
%      as: Md = M{1} + M{2}*xi(1) + M{3}*xi(2) + ... M{d+1}*xi(d)
%
% Additional notes:
%   - ALL the elements are assumed to have the same properties in terms
%     of MATERIAL and QUADRATURE rules.
%   ONLY if USEJULIA=1:
%   - this function uses the red_stiff_tensors_defects function 
%     implemented in the Julia module "DpROM.jl". 
%   - as such, this function supports ONLY models meshed with the elements
%     supported by both YetAnotherFEcode AND the DpROM.jl
%   - List of currently supported elements: 
%     Q4, Q8, TET4, TET10, HEX8, HEX20, WED15     (in YetAnotherFEcode)
%     Q4, Q8,       TET10, HEX8, HEX20, WED15     (in DpROM.jl)
%
% Reference: J. Marconi, P. Tiso, D.E. Quadrelli & F. Braghin, "A higher 
% order parametric nonlinear reduced order model for imperfect structures 
% using Neumann expansion", Nonlinear Dynamics, 2021.
%
% Created: 14 May 2021
% Last modified: 31 January 2022
% Author: Jacopo Marconi, Politecnico di Milano

function tensors = reduced_tensors_DpROM(myAssembly, elements, V, U, ...
    FORMULATION, volume, USEJULIA)

t0=tic;

% Handle incomplete input
if nargin < 7
    USEJULIA = 0;
elseif nargin < 6
    USEJULIA = 0;
    volume = 1;
elseif nargin < 5
    USEJULIA = 0;
    volume = 1;
    FORMULATION = 'N1';
end

% data from myAssembly
nodes    = myAssembly.Mesh.nodes;                   % nodes table
nel      = myAssembly.Mesh.nElements;           	% number of elements
nnodes   = myAssembly.Mesh.nNodes;               	% number of nodes
freedofs = myAssembly.Mesh.EBC.unconstrainedDOFs;   % free DOFs

% compute parametric mass matrix (always in matlab)
myMesh = myAssembly.Mesh;
DpROM = DpromAssembly(myMesh, U, V);
tensors.M = DpROM.ParametricMass;

if USEJULIA== 1
    % Element data (assuming ALL the elements have the same properties in 
    % terms of material and quadrature rule)
    Element = myAssembly.Mesh.Elements(1).Object;   % first element
    XGauss = Element.quadrature.X;                  % quadrature points
    WGauss = Element.quadrature.W;                  % quadrature weights
    C = Element.initialization.C;                   % constitutive matrix

    DIM = size(nodes,2);                            % 2D/3D problem
    ndofs = nnodes * DIM;                           % total number of DOFs
    nfree = length(freedofs);                       % number of free DOFs
    nd = size(U,2);                                 % number of defects

    % build connectivity table
    conn = zeros(nel, length(Element.iDOFs));
    for ee = 1 : nel
        conn(ee, :) = myAssembly.Mesh.Elements(ee).Object.iDOFs;
    end

    elements = int64(elements); % convert into integers
    conn = int64(conn);         % convert into integers

    % INPUT SIZE CHECK ____________________________________________________
    % for tensor assembly, always use the full DOFs displacement vectors
    if size(V,1)==nfree
        u = zeros(ndofs,size(V,2));
        u(freedofs,:) = V;
        V = u;
    elseif size(V,1)~=ndofs
        error('Wrong dimensions for V')
    end
    if size(U,1)==nfree
        u = zeros(ndofs,size(U,2));
        u(freedofs,:) = U;
        U = u;
    elseif isscalar(U) && U == 0
        U = V(:,1)*0;
    elseif size(U,1)~=ndofs
        error('Wrong dimensions for U')
    end
    FORMULATION = upper(FORMULATION);

    % JULIA _______________________________________________________________
    % add current path in Julia
    filename = 'wrapper\reduced_tensors_DpROM';         % this function name (up to folder "DpROM")
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

    % call the function once (for 1 element only) to precompile
    elem1 = elements(1,:);
    jl.call('red_stiff_tensors_defects', FORMULATION, volume, elem1, ...
        nodes, conn, C, V, U, XGauss, WGauss);

    disp([' REDUCED TENSORS (' FORMULATION ' ~ using Julia):'])
    fprintf(' Assembling %d elements ...', nel)

    % call the function in Julia to compute all the tensors
    a=jl.call('red_stiff_tensors_defects', FORMULATION, volume, ...
        elements, nodes, conn, C, V, U, XGauss, WGauss);

    % unpack results ______________________________________________________
    % depending on JULIA and juliamex versions the output may be different
    if iscell(a)
        Q2n= a{1}; %#ok<*GFLD>
        a2 = a{2};
        a3 = a{3};
        a4 = a{4};
        a5 = a{5};
        a6 = a{6};
        a7 = a{7};
        a8 = a{8};
        a9 = a{9};
        time = double(a{10})/1000;
    elseif isstruct(a)
        Q2n= getfield(a,'1'); %#ok<*GFLD>
        a2 = getfield(a,'2');
        a3 = getfield(a,'3');
        a4 = getfield(a,'4');
        a5 = getfield(a,'5');
        a6 = getfield(a,'6');
        a7 = getfield(a,'7');
        a8 = getfield(a,'8');
        a9 = getfield(a,'9');
        time = double(getfield(a,'10'))/1000;
    else
        error(' Output not recognized. Check Julia and juliamex versions')
    end
    if volume == 1
        ind = nd+1;
    else
        ind = 1; % (the other tensors are zero)
    end
    for ii = 1 : ind
        % if nd=1, squeeze singleton dimensions
        Q3d{ii}  = tensor(a2{ii}); %#ok<*AGROW>
        Q4dd{ii} = squeeze(tensor(a3{ii}));
        Q3n{ii}  = tensor(a4{ii});
        Q4d{ii}  = squeeze(permute(tensor(a5{ii}), [1 2 4 3]));    
        % if nd=1, DpROM.jl output already removes the final (5th) 
        % singleton dimension. Let us then distinguish:
        if nd>1
            Q5dd{ii} = permute(tensor(a6{ii}), [1 2 4 3 5]);
        else
            Q5dd{ii} = squeeze(permute(tensor(a6{ii}), [1 2 4 3]));
        end    
        Q4n{ii}  = tensor(a7{ii});
        Q5d{ii}  = squeeze(permute(tensor(a8{ii}), [1 2 3 5 4]));
        Q6dd{ii} = squeeze(permute(tensor(a9{ii}), [1 2 4 6 3 5]));
    end
        
    tensors.Q2n  = Q2n;
    tensors.Q3n  = Q3n;
    tensors.Q4n  = Q4n;
    tensors.Q3d  = Q3d;
    tensors.Q4d  = Q4d;
    tensors.Q5d  = Q5d;
    tensors.Q4dd = Q4dd;
    tensors.Q5dd = Q5dd;
    tensors.Q6dd = Q6dd;
    tensors.time = time;
    sftw = 'julia';
else
    disp([' REDUCED TENSORS (' FORMULATION ' ~ using Matlab):'])
    fprintf(' Assembling %d elements ...', nel)
    tensors = DpROM.Qtensors(FORMULATION, volume);
    tensors.M = DpROM.ParametricMass();
    time = toc(t0);
    sftw = 'matlab';
end

fprintf(' %.2f s (%.2f s)\n',toc(t0),time)
fprintf(' SPEED: %.1f el/s\n',nel/time)
fprintf(' SIZEs: %d - %d \n\n', size(V,2), size(U,2))

tensors.formulation = FORMULATION;
tensors.volume = volume;
tensors.software = sftw;

end


