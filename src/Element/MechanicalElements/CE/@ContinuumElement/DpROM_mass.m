% DpROM_mass
%
% CODE REFERENCE: ../YetAnotherFEcode/papers/DpROM
%
% Description: compute the defect-parametric reduced mass matrix using the
% approximated integral formulation used in [1] for the stiffness tensor.
%
% M: parametrized mass matrix, computed using an approximated
%    integration over the defected volume. The defected mass matrix can be
%    obtained as:
%           Md = M{1} + M{2}*xi(1) + M{3}*xi(2) + ... M{d+1}*xi(d)
%
% [1] Marconi et al. (2021). "A higher-order parametric nonlinear 
%     reduced-order model for imperfect structures using Neumann 
%     expansion". Nonlinear Dynamics. 
%     https://doi.org/10.1007/s11071-021-06496-y
%
% Created: 2021/10/21
% Jacopo Marconi, PhD, Politecnico di Milano

function M = DpROM_mass(self, Ve, Ue, Minit)
    M = Minit;
    nu = size(Ue,2);
    X = self.quadrature.X;
    W = self.quadrature.W;
    rho = self.Material.DENSITY;
    t = self.thickness;
    for ii = 1:length(W)
        Xi = X(:, ii);  % quadrature points
        we = W(ii);     % quadrature weights
        N = self.shape_functions(Xi);
        NN = kron(N', eye(self.nDim));
        [G, detJ] = shape_function_derivatives(self, Xi);
        Y = G*Ue;
        % integration of M through GAUSS QUADRATURE
        Mel_i = (NN'*NN)*(rho*t * we*detJ);
        Mel_i = Ve' * Mel_i * Ve;   % galerkin projection
        M{1}  = M{1} + Mel_i;
        % defect-additional mass matrices
        for d = 2 : nu + 1
            thd = Y(:, d-1);
            % approximated jacobian
            def_div = sum(thd(1 : self.nDim+1 : end));
            M{d}  = M{d}  + Mel_i  * def_div;
        end
    end
end