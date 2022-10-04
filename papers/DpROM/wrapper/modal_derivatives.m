% modal_derivatives
%
% Synthax:
% [MD, names] = modal_derivatives(myAssembly, elements, Phi, USEJULIA)
%
% Description: compute Modal Derivatives.
% INPUTS
%   - myAssembly: Assembly from YetAnotherFEcode. It MUST contain the field
%     ".DATA.K" storing the unconstrained linear stiffness matrix K0
%   - elements: table of the elements
%   - Phi: tall matrix containing by columns the selected Vibration Modes
%     (unconstrained)
%   - USEJULIA (optional): if set to 1, the Julia module "DpROM.jl" is used
%     (default value is 0).
% OUTPUTS
%   - MD: tall matrix containing by columns all the Modal Derivatives that
%     can be computed from the given Phi (unconstrained MDs are returned).
%   - names: matrix containing the subscripts of the MDs, for convenience.
%
% Additional notes:
%   - if USEJULIA=1, this function uses the stiffness_matrix_derivative 
%     function implemented in the Julia module "DpROM.jl". 
%   - if USEJULIA=1, this function supports ONLY models meshed with the 
%     elements supported by both YetAnotherFEcode AND the DpROM.jl
%   - List of currently supported elements: 
%     Q4, Q8, TET4, TET10, HEX8, HEX20, WED15   (in YetAnotherFEcode)
%     Q4, Q8,       TET10, HEX8, HEX20, WED15 	(in DpROM.jl)
%
% Created: 14 May 2021
% Author: Jacopo Marconi, Politecnico di Milano


function [MD, names] = modal_derivatives(myAssembly, elements, Phi, USEJULIA)

if nargin < 4
    USEJULIA = 0;
end

n = size(Phi,1);
n_VMs = size(Phi,2);

K0 = myAssembly.DATA.K;
K0 = myAssembly.constrain_matrix( K0 );

MD = zeros(n, n_VMs*(n_VMs+1)/2);
names = zeros(n_VMs*(n_VMs+1)/2, 2);
kk = 1;
for jj = 1 : n_VMs
    
    Phi_j = Phi(:, jj);
    if USEJULIA == 1
        dK_deta_j = stiffness_matrix_derivative(myAssembly, elements, Phi_j);
    else
        t0 = tic;
        fprintf(' dKdq, assembling %d elements ...', size(elements,1))
        dK_deta_j = myAssembly.matrix('stiffness_derivative',Phi_j);
        fprintf(' %.2f s\n',toc(t0))
    end
    dK_deta_j = myAssembly.constrain_matrix( dK_deta_j );
    
    for ii = 1 : n_VMs
        if ii < jj
            continue
        end
        
        Phi_i = myAssembly.constrain_vector( Phi(:, ii) );
        dPhi_i_deta_j = -K0\(dK_deta_j * Phi_i); 
        
        th =  dPhi_i_deta_j / max(abs(dPhi_i_deta_j));
        MD(:,kk) = myAssembly.unconstrain_vector( th );
        names(kk, :) = [ii jj];
        kk = kk + 1;
    end
end
