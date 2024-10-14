classdef Sensitivity
    % A class that describes the sensitivity of the Lyapunov Subcenter Manifolds (LSM).
    % Properties:
    %   sys: the system object.
    %   lsm: the Lyapunov Subcenter Manifolds object.
    %   domega0: sensitivity of the eigenfrequency.
    %   dphi0: sensitivity of the mode shape.
    %   LambdaTilde: normalized version of the state eigenvalue matrix.
    %   VTilde: normalized version of the state right eigenvector matrix.
    %   dVTilde: sensitivity of VTilde with respect to the eigenfrequency.
    %   UTilde: normalized version of the state left eigenvector matrix.
    %   dUTilde: sensitivity of UTilde with respect to the eigenfrequency.
    %   dR: cell matrix with the sensitivity of the R matrix (order x ndv).
    %   dW: cell matrix with the sensitivity of the W matrix (order x ndv).
    %   dgamma: cell matrix with the sensitivity of the gamma vector (order x ndv).
    % Methods:
    %   Sensitivity: constructor.
    %   sensitivity_backbone: sensitivity of a point on the backbone curve.
    %   sensitivity_lsm: sensitivity of the LSM.
    %   sensitivity_C: sensitivity of C.
    %   sensitivity_FoW: sensitivity of (F o W).
    %   compute_EN_tilde: compute E, NTilde, and dNTilde matrices.
    %   sensitivity_R: sensitivity of R.
    %   sensitivity_W: sensitivity of W.
    %   sensitivity_Y: sensitivity of Y.
    %   sensitivity_L: sensitivity of L.
    %   sensitivity_h: sensitivity of h.
    %   sensitivity_D: sensitivity of D.
    
    properties
        sys, lsm
        domega0, dphi0
        LambdaTilde, VTilde, dVTilde, UTilde, dUTilde
        dR, dW, dgamma
    end
    
    methods
        function obj = Sensitivity(sys, lsm)
            % Constructor.
            % Inputs:
            %   sys: system object.
            %   lsm: LSM object.

            % Store objects
            obj.sys = sys;
            obj.lsm = lsm;

            % Initialize system sensitivities
            obj.sys = obj.sys.second_order_form_sensitivity();
            obj.sys = obj.sys.state_space_sensitivity();

            % Initialize modal sensitivity
            [obj.dphi0, obj.domega0] = obj.sys.modal_analysis_sensitivity();

            % Initialize state space modal sensitivity
            [obj.LambdaTilde, ...
             obj.VTilde, obj.dVTilde, ...
             obj.UTilde, obj.dUTilde] = obj.sys.state_space_modal_analysis_sensitivity();

            % Initialize R, W, and gamma
            obj.dW = cell(lsm.order, obj.sys.ndv);
            obj.dR = cell(lsm.order, obj.sys.ndv);
            obj.dgamma = cell(lsm.order, obj.sys.ndv);

            % Order 1
            for ii = 1:obj.sys.ndv
                obj.dR{1, ii} = obj.domega0{ii} * obj.LambdaTilde;
                obj.dW{1, ii} = obj.domega0{ii} * kron(obj.dVTilde, obj.sys.phi0) ...
                                + kron(obj.VTilde, obj.dphi0{ii});
                obj.dgamma{1, ii} = obj.domega0{ii}; % * 1/2 * imag([1, -1] * sum(obj.dR{1, ii}, 2));
            end
        end

        function [Omega, dOmega] = sensitivity_backbone(obj, zTarget, iDof)
            % Sensitivity of a point on the backbone curve.
            % Inputs:
            %   zTarget: target point on the backbone curve.
            %   iDof: degree of freedom index.
            % Outputs:
            %   Omega: the frequency.
            %   dOmega: the sensitivity of the frequency.

            if zTarget == 0
                Omega = obj.sys.omega0;
                dOmega = [obj.dgamma{1, :}].';
            else
                % Find Rho and Omega
                Rho = obj.lsm.solve_rho(zTarget, iDof);
                Omega = obj.lsm.compute_omega(Rho);

                % Find zk
                [~, zk] = obj.lsm.compute_z_single(Rho, iDof);
    
                % Loop over theta
                num = zeros(obj.sys.ndv, 1);
                den = 0;
                for kk = 1:obj.lsm.nTheta
                    % Compute pk
                    pkTilde = exp(obj.lsm.theta(kk) .* [1j; -1j]);
    
                    % Loop over order
                    num_temp = zeros(obj.sys.ndv, 1);
                    den_temp = 0;
                    for ii = 1:obj.lsm.order
                        if ii == 1
                            pkTilde_ii = pkTilde;
                        else
                            pkTilde_ii = kron(pkTilde_ii, pkTilde);
                        end

                        % Loop over design variables
                        for jj = 1:obj.sys.ndv
                            num_temp(jj, 1) = num_temp(jj, 1) + Rho^ii * real(obj.dW{ii, jj}(iDof, :) * pkTilde_ii);
                        end
                        den_temp = den_temp + ii * Rho^(ii - 1) * real(obj.lsm.W{ii}(iDof, :) * pkTilde_ii);
                    end
    
                    % Update num and den
                    num = num + zk(kk) .* num_temp;
                    den = den + zk(kk) .* den_temp;
                end
    
                % Sensitivity of Rho
                dRho = -num./den;

                % Assemble sensitivity of Omega
                dOmega = 0;
                for ii = 2:obj.lsm.order
                    % Add terms that multiply the derivative of Rho
                    dOmega = dOmega + obj.lsm.gamma{ii} * (ii - 1) * Rho^(ii - 2);
                end
                dOmega = dOmega * dRho + [obj.dgamma{1, :}].';
                for ii = 2:obj.lsm.order
                    % Add terms related to the derivative of gamma
                    dOmega = dOmega + [obj.dgamma{ii, :}].' * Rho^(ii - 1);
                end
            end
        end

        function obj = sensitivity_lsm(obj)
            % Sensitivity of the Lyapunov Subcenter Manifolds (LSM).

            % Loop over orders
            for ii = 2:obj.lsm.order
                % Compute sensitivity of N
                [E, NTilde, dNTilde] = obj.compute_EN_tilde(ii);

                % Loop over design variables
                for jj = 1:obj.sys.ndv
                    % Sensitivity of C
                    dC = obj.sensitivity_C(ii, jj);
    
                    % Sensitivity of R
                    obj.dR{ii, jj} = obj.sensitivity_R(ii, jj, dC, E, NTilde, dNTilde);
    
                    % Sensitivity of gamma: 1/2 * imag([1, -1] * sum(obj.dR{ii, jj}, 2))
                    obj.dgamma{ii, jj} = full(imag(sum(obj.dR{ii, jj}(1, :))));
    
                    % Sensitivity of W
                    obj.dW{ii, jj} = obj.sensitivity_W(ii, jj, dC);
                end
            end
        end

        function dC = sensitivity_C(obj, order, dv)
            % Sensitivity of C
            % Inputs:
            %   order: current order.
            %   dv: design variable index.
            % Outputs:
            %   dC: the sensitivity of C.

            % Concatenate W and dW
            WdW = [obj.lsm.W, obj.dW(:, dv)];

            % Sensitivity of C
            dC = obj.sensitivity_FoW(order, WdW, dv);
            for jj = 2:(order - 1)
                [Y, dY] = obj.sensitivity_Y(order, jj, dv);
                dC = dC - obj.sys.dB{dv} * obj.lsm.W{jj} * Y ...
                        - obj.sys.B * obj.dW{jj, dv} * Y ...
                        - obj.sys.B * obj.lsm.W{jj} * dY;
            end
        end

        function dFoW = sensitivity_FoW(obj, order, WdW, dv)
            % Sensitivity of (F o W) at the given order.
            % Inputs:
            %   order: current order.
            %   WdW: cell array with the W and dW matrices involved in the computation.
            %   dv: design variable index.
            % Outputs:
            %   dFoW: the sensitivity of (F o W).
            
            % Derivative of F
            indices = obj.lsm.indicesFoW{order};
            dFoW = sparse(obj.sys.N, obj.lsm.Q^order);
            for ii = 1:length(indices) % loop over blocks
                idx = indices{ii};
                cols = size(idx, 2);
                if nnz(obj.sys.dF{cols, dv}) > 0
                    for jj = 1:size(idx, 1) % loop over indices
                        dFoW = dFoW + obj.sys.dF{cols, dv} * kronn(obj.lsm.W(idx(jj, :)));
                    end
                end
            end

            % Derivative of Wi
            for ii = 1:length(indices) 
                idx = indices{ii};
                cols = size(idx, 2);
                if nnz(obj.sys.F{cols}) > 0
                    for jj = 1:size(idx, 1)
                        for kk = 1:cols
                            idx_col = ones(1, cols);
                            idx_col(kk) = 2;
                            W_prod = 1;
                            for ll = 1:length(idx_col)
                                W_prod = kron(W_prod, WdW{idx(jj, ll), idx_col(ll)});
                            end
                            dFoW = dFoW + obj.sys.F{cols} * W_prod; % derivative of W
                        end
                    end
                end
            end
        end

        function [E, NTilde, dNTilde] = compute_EN_tilde(obj, order)
            % Compute E, NTilde, and dNTilde matrices of the normal for parameterization.
            % Inputs:
            %   order: current order.
            % Outputs:
            %   E: E matrix.
            %   NTilde: NTilde matrix.
            %   dNTilde: dNTilde matrix.

            % Find inner resonances
            [lIndices, jIndices] = find_inner_resonance(order);
            nInnerResonances = length(lIndices);

            % Initialize E and N triplets
            tripletsE = zeros(nInnerResonances, 3); % each column of E has nnz = 1
            tripletsN = zeros(nInnerResonances*2, 3); % each column of N has nnz = 2
            dNval = zeros(nInnerResonances*2, 1); % each column of N has nnz = 2

            % Loop over inner resonances
            itE = 0; itN = 0;
            for ii = 1:nInnerResonances

                % Triplets of E
                tripletsE(itE + 1, :) = [obj.lsm.Q*(lIndices(ii) - 1) + jIndices(ii), ii, 1];
                itE = itE + 1;

                % Triplets of NTilde
                tripletsN(itN + (1:2), :) = [2*(lIndices(ii) - 1) + (1:2).', ii*ones(2, 1), obj.UTilde(:, jIndices(ii))];

                % Triplets of dNTilde
                dNval(itN + (1:2), 1) = obj.dUTilde(:, jIndices(ii));
                itN = itN + 2;
            end

            % Build sparse matrices from triplets
            E = sparse(tripletsE(:, 1), tripletsE(:, 2), tripletsE(:, 3), obj.lsm.Q^order * obj.lsm.Q, nInnerResonances);
            NTilde = sparse(tripletsN(:, 1), tripletsN(:, 2), tripletsN(:, 3), obj.lsm.Q^order * 2, nInnerResonances);
            dNTilde = sparse(tripletsN(:, 1), tripletsN(:, 2), dNval, obj.lsm.Q^order * 2, nInnerResonances);
        end

        function dR = sensitivity_R(obj, order, dv, dC, E, NTilde, dNTilde)
            % Sensitivity of R
            % Inputs:
            %   order: current order.
            %   dv: design variable index.
            %   dC: sensitivity of C.
            %   E: E matrix.
            %   NTilde: NTilde matrix.
            %   dNTilde: dNTilde matrix.
            % Outputs:
            %   dR: the sensitivity of R.

            % Check order
            if mod(order, 2) == 0 % for even orders R = 0
                dR = sparse(obj.lsm.Q, obj.lsm.Q^order);
            else
                dN = obj.domega0{dv} * kron(dNTilde, obj.sys.phi0) ...
                     + kron(NTilde, obj.dphi0{dv});
                dRvec = E * (dN' * obj.lsm.C{order}(:) + kron(NTilde', obj.sys.phi0.') * dC(:));
                dR = reshape(dRvec, obj.lsm.Q, []);
            end
        end

        function dW = sensitivity_W(obj, order, dv, dC)
            % Sensitivity of W
            % Inputs:
            %   order: current order.
            %   dv: design variable index.
            %   dC: sensitivity of C.
            % Outputs:
            %   dW: the sensitivity of W.

            % Sensitivity of L
            [L, dL] = obj.sensitivity_L(order, dv);

            % Sensitivity of h
            [~, dh] = obj.sensitivity_h(order, dv, dC);

            % Sensitivity of W
            dWvec = lsqminnorm(L, dh - dL * obj.lsm.W{order}(:));
            dW = reshape(dWvec, obj.sys.N, []);
        end

        function [Y, dY] = sensitivity_Y(obj, ii, jj, dv)
            % Sensitivity of Y
            % Inputs:
            %   ii: row index (current order).
            %   jj: column index (previous orders).
            %   dv: design variable index.
            % Outputs:
            %   Y: the Y matrix.
            %   dY: the sensitivity of Y.

            Y = sparse(obj.lsm.Q^jj, obj.lsm.Q^ii);
            dY = sparse(obj.lsm.Q^jj, obj.lsm.Q^ii);
            thisR = obj.lsm.R{ii - jj + 1};
            thisdR = obj.dR{ii - jj + 1, dv};
            if nnz(thisR) > 0
                for kk = 1:jj
                    % Left identity
                    if kk > 1
                        leftI = speye(obj.lsm.Q^(kk - 1));
                    else
                        leftI = 1;
                    end
                    
                    % Right identity
                    if kk < jj
                        rightI = speye(obj.lsm.Q^(jj - kk));
                    else
                        rightI = 1;
                    end
    
                    % Build Y and dY
                    Y = Y + kron(kron(leftI, thisR), rightI); % kronn({leftI, thisR, rightI});
                    dY = dY + kron(kron(leftI, thisdR), rightI); % kronn({leftI, thisdR, rightI});
                end
            end
        end

        function [L, dL] = sensitivity_L(obj, order, dv)
            % Sensitivity of L
            % Inputs:
            %   order: current order.
            %   dv: design variable index.
            % Outputs:
            %   L: the L matrix.
            %   dL: the sensitivity of L.

            % Compute Y and dY
            [Y, dY] = obj.sensitivity_Y(order, order, dv);

            % Compute L and dL
            L = kron(Y.', obj.sys.B) - kron(speye(obj.lsm.Q^order), obj.sys.A);
            dL = kron(dY.', obj.sys.B) + kron(Y.', obj.sys.dB{dv}) ...
                 - kron(speye(obj.lsm.Q^order), obj.sys.dA{dv});
        end

        function [h, dh] = sensitivity_h(obj, order, dv, dC)
            % Sensitivity of h
            % Inputs:
            %   order: current order.
            %   dv: design variable index.
            %   dC: sensitivity of C.
            % Outputs:
            %   h: the h vector.
            %   dh: the sensitivity of h.

            % Compute D and dD
            [D, dD] = obj.sensitivity_D(order, dv);

            % Compute h and dh
            h = obj.lsm.C{order}(:) - D * obj.lsm.R{order}(:);
            dh = dC(:) - dD * obj.lsm.R{order}(:) - D * obj.dR{order, dv}(:);
        end

        function [D, dD] = sensitivity_D(obj, order, dv)
            % Sensitivity of D
            % Inputs:
            %   order: current order.
            %   dv: design variable index.
            % Outputs:
            %   D: the D matrix.
            %   dD: the sensitivity of D.

            D = kron(speye(obj.lsm.Q^order), obj.sys.B*obj.lsm.W{1});
            dD = kron(speye(obj.lsm.Q^order), obj.sys.dB{dv}*obj.lsm.W{1} + obj.sys.B*obj.dW{1, dv});
        end
    end
end
