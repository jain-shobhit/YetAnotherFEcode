% EXAMPLE: beam meshed with multiple (yet compatible) element types
clear;
close all;
clc

%% PREPARE MODEL                                                    

% Material
E       = 70e9;     % Young's modulus [Pa]
rho     = 2700;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio

myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);

% MESH_____________________________________________________________________
filename = 'Job-BeamMixed';
meshinfo = abqmesh(filename);
nodes    = meshinfo.nodes;
elements = meshinfo.elem;
nset     = meshinfo.nset;

% Element constructor (define one for each element type)
myElementConstructor{1} = @()Wed15Element(myMaterial);
myElementConstructor{2} = @()Tet10Element(myMaterial);

% Build Mesh object
myMesh = Mesh(nodes);
myMesh.create_elements_table(elements, myElementConstructor);

% MESH > BOUNDARY CONDITONS
myMesh.set_essential_boundary_condition([nset{1} nset{2}],1:3,0)

% ASSEMBLY ________________________________________________________________
myAssembly = Assembly(myMesh);
M = myAssembly.mass_matrix();
nNodes = size(nodes,1);
u0 = zeros( myMesh.nDOFs, 1);
[K,~] = myAssembly.tangent_stiffness_and_force(u0);

% store matrices
myAssembly.DATA.K = K;
myAssembly.DATA.M = M;


%% EXAMPLE 1: vibration modes                                       

% Eigenvalue problem_______________________________________________________
n_VMs = 3; % first n_VMs modes with lowest frequency calculated
Kc = myAssembly.constrain_matrix(K);
Mc = myAssembly.constrain_matrix(M);
[V0,om] = eigs(Kc,Mc,n_VMs,'SM');
[f0,ind] = sort(sqrt(diag(om))/2/pi);
V0 = V0(:,ind);
for ii = 1:n_VMs
    V0(:,ii) = V0(:,ii)/max(sqrt(sum(V0(:,ii).^2,2)));
end
V0 = myAssembly.unconstrain_vector(V0);

% PLOT
mod = 1;
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMesh(nodes,elements,0);
v1 = reshape(V0(:,mod),3,[]).';
PlotFieldonDeformedMesh(nodes,elements,v1,'factor',10);
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0(mod),3) ' Hz'])
drawnow

%% EXAMPLE 2: static test                                           

% Define external force:
% Body force
Pressure = 1e6;
F = Pressure * myAssembly.uniform_body_force();

[u_lin, u] = static_equilibrium(myAssembly, F, 'method', 'newton', ...
    'nsteps', 10, 'display', 'iter');
ULIN = reshape(u_lin,3,[]).';	% Linear response
UNL = reshape(u,3,[]).';        % Nonlinear response

[K,f] = myAssembly.tangent_stiffness_and_force(u);

fprintf(['\n <strong>Max displacements</strong>:\n  Linear:\t\t%.3i \n' ...
    '  Nonlinear:\t%.3i \n\n'],max(u_lin(:)),max(u(:)))

% PLOT
figure('units','normalized','position',[.2 .1 .6 .8])
scale = 1;
PlotMesh(nodes,elements,0);
PlotFieldonDeformedMesh(nodes,elements,UNL,'factor',scale,'color','w');
PlotFieldonDeformedMesh(nodes,elements,ULIN,'factor',scale,'color','w');
colormap jet
title(['NONLINEAR STATIC RESPONSE (scale factor: ' num2str(scale) 'x)'])
axis on
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]'); grid on
legend('mesh','NL','LIN','location','northeast')

%% EXAMPLE 3: compute modal derivatives                             

% take the first eigenmode
Phi1 = V0(:,1);
Phi1_c = myAssembly.constrain_vector(Phi1);

% compute dK/dq1
dK = myAssembly.stiffness_derivative(Phi1);
dK_c = myAssembly.constrain_matrix(dK);

% compute MD11
MD11_c = -Kc\(dK_c * Phi1_c);
MD11 = myAssembly.unconstrain_vector(MD11_c);
MD11 = -MD11/max(abs(MD11(:))); % normalize to 1

% PLOT
figure('units','normalized','position',[.2 .1 .6 .8])
v1 = reshape(MD11,3,[]).';
PlotFieldonDeformedMesh(nodes,elements,v1,'factor',2,'component','U1');
title('MD_{11} (U1)')
drawnow


%% EXAMPLE 4: reduced stiffness tensors                             

RB = [V0(:,1) MD11]; % reduced basis with VM_1 and MD_11
m = size(RB,2);

RmyAssembly = ReducedAssembly(myMesh, RB);

% reduced stiffness tensors
K2r = RB'*K*RB;
ndofs_per_element = myMesh.Elements(1).Object.nNodes * myMesh.Elements(1).Object.nDim;
if m > ndofs_per_element 
    % compute element tensor (size ndofs_per_element), then project
    mode = 'standard';
    disp(' Standard tensor construction:')
    
    tic
    K3r = RmyAssembly.tensor('T2',[m m m], [2 3], mode);
    fprintf(' K3r... %.2f s (%.1f elem/s)\n',toc,myMesh.nElements/toc)
    
    tic
    K4r = RmyAssembly.tensor('T3',[m m m m], [2 3 4], mode);
    fprintf(' K4r... %.2f s (%.1f elem/s)\n',toc,myMesh.nElements/toc)
else
    % use Element-Level Projection (directly computes the projected tensor, of size m)
    mode = 'ELP';
    disp(' Element-wise projection:')
    
    tic
    K3r = RmyAssembly.tensor('T2',[m m m], [2 3], mode);
    fprintf(' K3r... %.2f s (%.1f elem/s)\n',toc,myMesh.nElements/toc)

    tic
    K4r = RmyAssembly.tensor('T3',[m m m m], [2 3 4], mode);
    fprintf(' K4r... %.2f s (%.1f elem/s)\n',toc,myMesh.nElements/toc)
end

% compute tensors for the tangent stiffness matrix
K3rt = K3r + permute(K3r, [1 3 2]); 
K4rt = K4r + permute(K4r, [1 3 2 4]) + permute(K4r, [1 4 2 3]);



