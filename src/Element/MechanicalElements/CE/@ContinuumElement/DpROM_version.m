% DpROM_version
%
% CODE REFERENCE: ../YetAnotherFEcode/papers/DpROM
%
% Reference: 
%   Marconi et al. (2021). "A higher-order parametric nonlinear 
%   reduced-order model for imperfect structures using Neumann 
%   expansion". Nonlinear Dynamics. 
%   https://doi.org/10.1007/s11071-021-06496-y
%
% DpROM versions:
%   -N0n:  0-th order Neumann expansion, integration over nominal volume
%   -N0d:  0-th order Neumann expansion, integration over defected volume
%   -N1Tn: 1-th order Neumann expansion, integration over nominal volume*
%   -N1Td: 1-th order Neumann expansion, integration over defected volume*
%   -N1n:  1-th order Neumann expansion, integration over nominal volume
%   -N1d:  1-th order Neumann expansion, integration over defected volume
% *truncated version (no 5-th and 6-th order tensors)
%
% Created: 2021/12/13
% Jacopo Marconi, PhD, Politecnico di Milano

function Q = DpROM_version(self, Ve, Ue, data, fun_name)
    Q = feval(fun_name, self, Ve, Ue, data);    
end

function Q = Qten_N0n(self, Ve, Ue, data)
% DpROM-N0n: zeroth order Neumann expansion, integration over nominal volume

    L1 = data.L1;
    C0 = self.initialization.C;
    H = self.initialization.H;
    Q = data.Qinit;

    X = self.quadrature.X;
    W = self.quadrature.W;
    for ii = 1:length(W)
        Xi = X(:, ii);  % quadrature points
        we = W(ii);     % quadrature weights
        [G, detJ] = shape_function_derivatives(self, Xi);               
        Y = G*Ue;
        G = G*Ve;
        C = C0*(we * detJ); % <--- includes we and detJ!

        % Q2n, nominal
        HG = H*G;
        Q2n_i = HG'*C*HG;

        % Q3dd, defected
        CHG = C*HG;
        L1YG = einsum('kla,aK,lJ->kJK',L1,Y,G);
        Q3d_i = einsum('kI,kJK->IJK',CHG,L1YG);
        Q3d_i = Q3d_i + permute(Q3d_i,[2 1 3]);

        % Q4dd, defected
        Q4dd_i = einsum('jIK,jk,kJL->IJKL',L1YG,C,L1YG);

        % Q3n, nominal 
        L1GG = einsum('kla,aK,lJ->kJK',L1,G,G);
        Q3n_i = einsum('kI,kJK->IJK',CHG,L1GG);
        Q3n_i = 0.5 * Q3n_i + permute(Q3n_i,[2 1 3]);

        % Q4d, defected
        Q4d_a = einsum('jIL,jk,kJK->IJKL',L1YG,C,L1GG);
        Q4d_i = 0.5*Q4d_a + permute(Q4d_a,[3 2 1 4]);

        % Q4n, nominal
        Q4n_i = 0.5 * einsum('jIK,jk,kLJ->IJKL',L1GG,C,L1GG);

        Q.Q2n{1}  = Q.Q2n{1}  + Q2n_i;
        Q.Q3d{1}  = Q.Q3d{1}  + Q3d_i;
        Q.Q4dd{1} = Q.Q4dd{1} + Q4dd_i;
        Q.Q3n{1}  = Q.Q3n{1}  + Q3n_i;
        Q.Q4d{1}  = Q.Q4d{1}  + Q4d_i;
        Q.Q4n{1}  = Q.Q4n{1}  + Q4n_i;
    end
end

function Q = Qten_N0d(self, Ve, Ue, data)
% DpROM-N0d: zeroth order Neumann expansion, integration over defected volume

    L1 = data.L1;
    C0 = self.initialization.C;
    H = self.initialization.H;
    Q = data.Qinit;
    nu = size(Ue,2);

    X = self.quadrature.X;
    W = self.quadrature.W;
    for ii = 1:length(W)
        Xi = X(:, ii);  % quadrature points
        we = W(ii);     % quadrature weights
        [G, detJ] = shape_function_derivatives(self, Xi);               
        Y = G*Ue;
        G = G*Ve;
        C = C0*(we * detJ); % <--- includes we and detJ!

        % Q2n, nominal
        HG = H*G;
        Q2n_i = HG'*C*HG;

        % Q3dd, defected
        CHG = C*HG;
        L1YG = einsum('kla,aK,lJ->kJK',L1,Y,G);
        Q3d_i = einsum('kI,kJK->IJK',CHG,L1YG);
        Q3d_i = Q3d_i + permute(Q3d_i,[2 1 3]);

        % Q4dd, defected
        Q4dd_i = einsum('jIK,jk,kJL->IJKL',L1YG,C,L1YG);

        % Q3n, nominal 
        L1GG = einsum('kla,aK,lJ->kJK',L1,G,G);
        Q3n_i = einsum('kI,kJK->IJK',CHG,L1GG);
        Q3n_i = 0.5 * Q3n_i + permute(Q3n_i,[2 1 3]);

        % Q4d, defected
        Q4d_a = einsum('jIL,jk,kJK->IJKL',L1YG,C,L1GG);
        Q4d_i = 0.5*Q4d_a + permute(Q4d_a,[3 2 1 4]);

        % Q4n, nominal
        Q4n_i = 0.5 * einsum('jIK,jk,kLJ->IJKL',L1GG,C,L1GG);

        Q.Q2n{1}  = Q.Q2n{1}  + Q2n_i;
        Q.Q3d{1}  = Q.Q3d{1}  + Q3d_i;
        Q.Q4dd{1} = Q.Q4dd{1} + Q4dd_i;
        Q.Q3n{1}  = Q.Q3n{1}  + Q3n_i;
        Q.Q4d{1}  = Q.Q4d{1}  + Q4d_i;
        Q.Q4n{1}  = Q.Q4n{1}  + Q4n_i;

        for d = 2 : nu + 1
            thd = Y(:, d-1);
            def_div = sum(thd(1 : self.nDim+1 : end));
            Q.Q2n{d}  = Q.Q2n{d}  + Q2n_i  * def_div;
            Q.Q3d{d}  = Q.Q3d{d}  + Q3d_i  * def_div;
            Q.Q4dd{d} = Q.Q4dd{d} + Q4dd_i * def_div;
            Q.Q3n{d}  = Q.Q3n{d}  + Q3n_i  * def_div;
            Q.Q4d{d}  = Q.Q4d{d}  + Q4d_i  * def_div;
            Q.Q4n{d}  = Q.Q4n{d}  + Q4n_i  * def_div;
        end
    end
end

function Q = Qten_N1Tn(self, Ve, Ue, data)
% DpROM-N1Tn: first order Neumann expansion, integration over nominal 
%             volume (truncated version - no fifth and sixth order tensors)

    L1 = data.L1;
    L2 = data.L2;          
    C0 = self.initialization.C;
    H = self.initialization.H;
    Q = data.Qinit;

    X = self.quadrature.X;
    W = self.quadrature.W;
    for ii = 1:length(W)
        Xi = X(:, ii);  % quadrature points
        we = W(ii);     % quadrature weights
        [G, detJ] = shape_function_derivatives(self, Xi);               
        Y = G*Ue;
        G = G*Ve;
        C = C0*(we * detJ); % <--- includes we and detJ!

        % Q2n, nominal
        HG = H*G;
        Q2n_i = HG'*C*HG;

        % Q3dd, defected
        CHG = C*HG;
        L2YG = einsum('kla,aK,lJ->kJK',L2,Y,G);
        Q3d_i = einsum('kI,kJK->IJK',CHG,L2YG);
        Q3d_i = Q3d_i + permute(Q3d_i,[2 1 3]);

        % Q4dd, defected
        Q4dd_i = einsum('jIK,jk,kJL->IJKL',L2YG,C,L2YG);

        % Q3n, nominal 
        L1GG = einsum('kla,aK,lJ->kJK',L1,G,G);
        Q3n_i = einsum('kI,kJK->IJK',CHG,L1GG);
        Q3n_i = 0.5 * Q3n_i + permute(Q3n_i,[2 1 3]);

        % Q4d, defected
        Q4d_a = einsum('jIL,jk,kJK->IJKL',L2YG,C,L1GG);
        Q4d_i = 0.5*Q4d_a + permute(Q4d_a,[3 2 1 4]);

        % Q4n, nominal
        Q4n_i = 0.5 * einsum('jIK,jk,kLJ->IJKL',L1GG,C,L1GG);

        Q.Q2n{1}  = Q.Q2n{1}  + Q2n_i;
        Q.Q3d{1}  = Q.Q3d{1}  + Q3d_i;
        Q.Q4dd{1} = Q.Q4dd{1} + Q4dd_i;
        Q.Q3n{1}  = Q.Q3n{1}  + Q3n_i;
        Q.Q4d{1}  = Q.Q4d{1}  + Q4d_i;
        Q.Q4n{1}  = Q.Q4n{1}  + Q4n_i;
    end
end

function Q = Qten_N1Td(self, Ve, Ue, data)
% DpROM-N1Td: first order Neumann expansion, integration over defected 
%             volume (truncated version - no fifth and sixth order tensors)

    L1 = data.L1;
    L2 = data.L2;           
    C0 = self.initialization.C;
    H = self.initialization.H;
    Q = data.Qinit;
    nu = size(Ue,2);

    X = self.quadrature.X;
    W = self.quadrature.W;
    for ii = 1:length(W)
        Xi = X(:, ii);  % quadrature points
        we = W(ii);     % quadrature weights
        [G, detJ] = shape_function_derivatives(self, Xi);               
        Y = G*Ue;
        G = G*Ve;
        C = C0*(we * detJ); % <--- includes we and detJ!

        % Q2n, nominal
        HG = H*G;
        Q2n_i = HG'*C*HG;

        % Q3dd, defected
        CHG = C*HG;
        L2YG = einsum('kla,aK,lJ->kJK',L2,Y,G);
        Q3d_i = einsum('kI,kJK->IJK',CHG,L2YG);
        Q3d_i = Q3d_i + permute(Q3d_i,[2 1 3]);

        % Q4dd, defected
        Q4dd_i = einsum('jIK,jk,kJL->IJKL',L2YG,C,L2YG);

        % Q3n, nominal 
        L1GG = einsum('kla,aK,lJ->kJK',L1,G,G);
        Q3n_i = einsum('kI,kJK->IJK',CHG,L1GG);
        Q3n_i = 0.5 * Q3n_i + permute(Q3n_i,[2 1 3]);

        % Q4d, defected
        Q4d_a = einsum('jIL,jk,kJK->IJKL',L2YG,C,L1GG);
        Q4d_i = 0.5*Q4d_a + permute(Q4d_a,[3 2 1 4]);

        % Q4n, nominal
        Q4n_i = 0.5 * einsum('jIK,jk,kLJ->IJKL',L1GG,C,L1GG);

        Q.Q2n{1}  = Q.Q2n{1}  + Q2n_i;
        Q.Q3d{1}  = Q.Q3d{1}  + Q3d_i;
        Q.Q4dd{1} = Q.Q4dd{1} + Q4dd_i;
        Q.Q3n{1}  = Q.Q3n{1}  + Q3n_i;
        Q.Q4d{1}  = Q.Q4d{1}  + Q4d_i;
        Q.Q4n{1}  = Q.Q4n{1}  + Q4n_i;

        for d = 2 : nu + 1
            thd = Y(:, d-1);
            def_div = sum(thd(1 : self.nDim+1 : end));
            Q.Q2n{d}  = Q.Q2n{d}  + Q2n_i  * def_div;
            Q.Q3d{d}  = Q.Q3d{d}  + Q3d_i  * def_div;
            Q.Q4dd{d} = Q.Q4dd{d} + Q4dd_i * def_div;
            Q.Q3n{d}  = Q.Q3n{d}  + Q3n_i  * def_div;
            Q.Q4d{d}  = Q.Q4d{d}  + Q4d_i  * def_div;
            Q.Q4n{d}  = Q.Q4n{d}  + Q4n_i  * def_div;
        end
    end
end

function Q = Qten_N1n(self, Ve, Ue, data)
% DpROM-N1n: first order Neumann expansion, integration over nominal volume

    L1 = data.L1;
    L2 = data.L2;
    L3 = data.L3;            
    C0 = self.initialization.C;
    H = self.initialization.H;
    Q = data.Qinit;
    nu = size(Ue,2);

    X = self.quadrature.X;
    W = self.quadrature.W;
    for ii = 1:length(W)
        Xi = X(:, ii);  % quadrature points
        we = W(ii);     % quadrature weights
        [G, detJ] = shape_function_derivatives(self, Xi);               
        Y = G*Ue;
        G = G*Ve;
        C = C0*(we * detJ); % <--- includes we and detJ!

        % Q2n, nominal
        HG = H*G;
        Q2n_i = HG'*C*HG;

        % Q3dd, defected
        CHG = C*HG;
        L2YG = einsum('kla,aK,lJ->kJK',L2,Y,G);
        Q3d_i = einsum('kI,kJK->IJK',CHG,L2YG);
        Q3d_i = Q3d_i + permute(Q3d_i,[2 1 3]);

        % Q4dd, defected
        Q4dd_i = einsum('jIK,jk,kJL->IJKL',L2YG,C,L2YG);

        % Q3n, nominal 
        L1GG = einsum('kla,aK,lJ->kJK',L1,G,G);
        Q3n_i = einsum('kI,kJK->IJK',CHG,L1GG);
        Q3n_i = 0.5 * Q3n_i + permute(Q3n_i,[2 1 3]);

        % Q4d, defected
        Q4d_a = einsum('jIL,jk,kJK->IJKL',L2YG,C,L1GG);
        Q4d_a = 0.5*Q4d_a + permute(Q4d_a,[3 2 1 4]);
        L3YGG = einsum('klab,bL,aK,lJ->kJKL',L3,Y,G,G);
        Q4d_b = einsum('kI,kJKL->IJKL',CHG,L3YGG);
        Q4d_b = 2*permute(Q4d_b, [3 2 1 4]) + Q4d_b;
        Q4d_i = Q4d_a + Q4d_b;

        % Q5dd, defected
        Q5dd_i = einsum('jIL,jk,kJKM->IJKLM',L2YG,C,L3YGG);
        Q5dd_i = 2*permute(Q5dd_i, [3 2 1 5 4]) + Q5dd_i;

        % Q4n, nominal
        Q4n_i = 0.5 * einsum('jIK,jk,kLJ->IJKL',L1GG,C,L1GG);

        % Q5d, defected
        Q5d_i = einsum('jIK,jk,kJLM->IKJLM',L1GG,C,L3YGG);
        Q5d_i = Q5d_i + permute(Q5d_i, [4 3 2 1 5]);

        % Q6dd, defected
        Q6dd_i = 2*einsum('jILM,jk,kJKN->IKJLMN',L3YGG,C,L3YGG);

        Q.Q2n{1}  = Q.Q2n{1}  + Q2n_i;
        Q.Q3d{1}  = Q.Q3d{1}  + Q3d_i;
        Q.Q4dd{1} = Q.Q4dd{1} + Q4dd_i;
        Q.Q3n{1}  = Q.Q3n{1}  + Q3n_i;
        Q.Q4d{1}  = Q.Q4d{1}  + Q4d_i;
        Q.Q5dd{1} = Q.Q5dd{1} + Q5dd_i;
        Q.Q4n{1}  = Q.Q4n{1}  + Q4n_i;
        Q.Q5d{1}  = Q.Q5d{1}  + Q5d_i;
        Q.Q6dd{1} = Q.Q6dd{1} + Q6dd_i;
    end
end

function Q = Qten_N1d(self, Ve, Ue, data)
% DpROM-N1n: first order Neumann expansion, integration over defected volume

    L1 = data.L1;
    L2 = data.L2;
    L3 = data.L3;            
    C0 = self.initialization.C;
    H = self.initialization.H;
    Q = data.Qinit;
    nu = size(Ue,2);

    X = self.quadrature.X;
    W = self.quadrature.W;
    for ii = 1:length(W)
        Xi = X(:, ii);  % quadrature points
        we = W(ii);     % quadrature weights
        [G, detJ] = shape_function_derivatives(self, Xi);               
        Y = G*Ue;
        G = G*Ve;
        C = C0*(we * detJ); % <--- includes we and detJ!

        % Q2n, nominal
        HG = H*G;
        Q2n_i = HG'*C*HG;

        % Q3dd, defected
        CHG = C*HG;
        L2YG = einsum('kla,aK,lJ->kJK',L2,Y,G);
        Q3d_i = einsum('kI,kJK->IJK',CHG,L2YG);
        Q3d_i = Q3d_i + permute(Q3d_i,[2 1 3]);

        % Q4dd, defected
        Q4dd_i = einsum('jIK,jk,kJL->IJKL',L2YG,C,L2YG);

        % Q3n, nominal 
        L1GG = einsum('kla,aK,lJ->kJK',L1,G,G);
        Q3n_i = einsum('kI,kJK->IJK',CHG,L1GG);
        Q3n_i = 0.5 * Q3n_i + permute(Q3n_i,[2 1 3]);

        % Q4d, defected
        Q4d_a = einsum('jIL,jk,kJK->IJKL',L2YG,C,L1GG);
        Q4d_a = 0.5*Q4d_a + permute(Q4d_a,[3 2 1 4]);
        L3YGG = einsum('klab,bL,aK,lJ->kJKL',L3,Y,G,G);
        Q4d_b = einsum('kI,kJKL->IJKL',CHG,L3YGG);
        Q4d_b = 2*permute(Q4d_b, [3 2 1 4]) + Q4d_b;
        Q4d_i = Q4d_a + Q4d_b;

        % Q5dd, defected
        Q5dd_i = einsum('jIL,jk,kJKM->IJKLM',L2YG,C,L3YGG);
        Q5dd_i = 2*permute(Q5dd_i, [3 2 1 5 4]) + Q5dd_i;

        % Q4n, nominal
        Q4n_i = 0.5 * einsum('jIK,jk,kLJ->IJKL',L1GG,C,L1GG);

        % Q5d, defected
        Q5d_i = einsum('jIK,jk,kJLM->IKJLM',L1GG,C,L3YGG);
        Q5d_i = Q5d_i + permute(Q5d_i, [4 3 2 1 5]);

        % Q6dd, defected
        Q6dd_i = 2*einsum('jILM,jk,kJKN->IKJLMN',L3YGG,C,L3YGG);

        Q.Q2n{1}  = Q.Q2n{1}  + Q2n_i;
        Q.Q3d{1}  = Q.Q3d{1}  + Q3d_i;
        Q.Q4dd{1} = Q.Q4dd{1} + Q4dd_i;
        Q.Q3n{1}  = Q.Q3n{1}  + Q3n_i;
        Q.Q4d{1}  = Q.Q4d{1}  + Q4d_i;
        Q.Q5dd{1} = Q.Q5dd{1} + Q5dd_i;
        Q.Q4n{1}  = Q.Q4n{1}  + Q4n_i;
        Q.Q5d{1}  = Q.Q5d{1}  + Q5d_i;
        Q.Q6dd{1} = Q.Q6dd{1} + Q6dd_i;

        for d = 2 : nu + 1
            thd = Y(:, d-1);
            def_div = sum(thd(1 : self.nDim+1 : end));
            Q.Q2n{d}  = Q.Q2n{d}  + Q2n_i  * def_div;
            Q.Q3d{d}  = Q.Q3d{d}  + Q3d_i  * def_div;
            Q.Q4dd{d} = Q.Q4dd{d} + Q4dd_i * def_div;
            Q.Q3n{d}  = Q.Q3n{d}  + Q3n_i  * def_div;
            Q.Q4d{d}  = Q.Q4d{d}  + Q4d_i  * def_div;
            Q.Q5dd{d} = Q.Q5dd{d} + Q5dd_i * def_div;
            Q.Q4n{d}  = Q.Q4n{d}  + Q4n_i  * def_div;
            Q.Q5d{d}  = Q.Q5d{d}  + Q5d_i  * def_div;
            Q.Q6dd{d} = Q.Q6dd{d} + Q6dd_i * def_div;
        end
    end
end