% defect_sensitivities
%
% Synthax:
% [DS, names] = defect_sensitivities(myAssembly, elements, Phi, U, formulation, USEJULIA)
%
% Description: compute Defect Sensitivities.
% INPUTS
%   - myAssembly: Assembly from YetAnotherFEcode. It MUST contain the field
%     ".DATA.K" storing the unconstrained linear stiffness matrix K0
%   - elements: table of the elements
%   - Phi: tall matrix collecting by columns the selected Vibration Modes
%     (unconstrained)
%   - U: tall matrix collecting by columns the (unconstrained) defect shapes
%   - formulation (optional): choose the Neumann formulation between "N1",
%     "N1T" and "N0" (N1 and N1T are the same here).
%     The default value is N1.
%   - USEJULIA (optional): if set to 1, the Julia module "DpROM.jl" is used
%     (default value is 0).
% OUTPUTS
%   - DS: tall matrix containing by columns all the Defect Sensitivities
%     that can be computed from the given Phi and U (unconstrained DSs are 
%     returned).
%   - names: matrix containing the subscripts of the DSs, for convenience.
%
% Additional notes:
%   - if USEJULIA=1, this function uses the stiffness_matrix_sensitivity 
%     function implemented in the Julia module "DpROM.jl". 
%   - if USEJULIA=1, this function supports ONLY models meshed with the 
%     elements supported by both YetAnotherFEcode AND the DpROM.jl
%   - List of currently supported elements: 
%     Q4, Q8, TET4, TET10, HEX8, HEX20, WED15   (in YetAnotherFEcode)
%     Q4, Q8, TET4, TET10, HEX8, HEX20, WED15 	(in DpROM.jl)
%
% Reference: J. Marconi, P. Tiso, D.E. Quadrelli & F. Braghin, "A higher 
% order parametric nonlinear reduced order model for imperfect structures 
% using Neumann expansion", Nonlinear Dynamics, 2021.
%
% Created: 14 May 2021
% Author: Jacopo Marconi, Politecnico di Milano

function [DS, names] = defect_sensitivities(myAssembly, elements, Phi, U, formulation, USEJULIA)

if nargin < 5
    formulation = 'N1';
    USEJULIA = 0;
    fprintf(' Default formulation for DS is: N1')
elseif nargin < 6
    USEJULIA = 0;
end

n  = size(Phi, 1);
nm = size(Phi, 2);
nd = size(U,   2);

K0 = myAssembly.DATA.K;
K0 = myAssembly.constrain_matrix( K0 );

DS = zeros(n, nm*nd);
names = zeros(nm*nd, 2);
cc = 1;
for kk = 1 : nd
    Uk = U(:,kk);  
    if USEJULIA == 1
        dK_dxi_k = stiffness_matrix_sensitivity(myAssembly, elements, Uk, formulation);
    else
        t0 = tic;
        fprintf(' dKdxi, assembling %d elements ...', size(elements,1))
        dK_dxi_k = myAssembly.matrix('DpROM_stiffness_defect_derivative',Uk, formulation);
        fprintf(' %.2f s\n',toc(t0))
    end
    dK_dxi_k = myAssembly.constrain_matrix( dK_dxi_k );
    for ii = 1 : nm
        Phi_i = myAssembly.constrain_vector( Phi(:, ii) );
        Xi = - K0 \ (dK_dxi_k * Phi_i);
        DS(:, cc) =  myAssembly.unconstrain_vector( Xi ) / max(abs(Xi));
        names(cc, :) = [ii kk]';
        cc = cc+1;
    end
end
disp(' ')