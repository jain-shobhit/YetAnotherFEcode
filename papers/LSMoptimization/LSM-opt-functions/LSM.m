classdef LSM
    % A class that handles the Lyapunov Subcenter Manifolds (LSM).
    % Properties:
    %   sys: system object.
    %   order: order of the LSM.
    %   Q: number of modes in the master subspace.
    %   W: cell array with the coefficients of the manifold parametrization.
    %   R: cell array with the coefficients of the normal form parametrization.
    %   gamma: vector of the gamma values.
    %   C: cell array with the coefficients of the nonlinearity.
    %   nTheta: number of theta values.
    %   theta: vector of theta values.
    %   indicesFoW: cell array with the indices for the (F o W) computation.
    % Methods:
    %   compute_manifold: compute the LSM manifold.
    %   compute_Y: compute the Y matrix for the computation of the C matrix.
    %   compute_C: compute the C matrix for the computation of the W matrix.
    %   get_indices_FoW: get the indices for the computation of (F o W).
    %   compute_FoW: compute the (F o W) matrix.
    %   compute_EN: compute the E and N matrices for the computation of the R matrix.
    %   normal_form_parametrization: compute the normal form parametrization coefficients.
    %   compute_L: compute the L matrix.
    %   compute_D: compute the D matrix.
    %   compute_h: compute the h vector.
    %   solve_cohomological: solve the cohomological equation.
    %   compute_omega: compute the omega values.
    %   compute_z_single: compute the rms value of z for a single rho value.
    %   compute_z: compute the rms value of z for a range of rho values.
    %   nonlcon: nonlinear equality constraints for fmincon.
    %   solve_rho_at_leading_order: solve for the rho value at the leading order.
    %   solve_rho: solve for the rho values given the rms values of z.

    properties
        sys
        order, Q
        W, R, gamma
        C
        nTheta, theta
        indicesFoW
    end
    
    methods
        function obj = LSM(sys, order, nTheta)
            % Constructor.
            % Inputs:
            %   sys: system object.
            %   order: order of the LSM.
            %   nTheta: number of theta values. Default is 2^5.

            % System object
            obj.sys = sys;

            % Check number of inputs
            if nargin == 2
                obj.nTheta = 2^5;
            else
                obj.nTheta = nTheta;
            end
            obj.theta = linspace(0, 2*pi, obj.nTheta);

            % Max expansion order
            obj.Q = 2;
            obj.order = order;

            % Initialize variables
            obj.C = cell(order, 1);
            obj.R = cell(order, 1);
            obj.gamma = cell(order, 1);
            obj.W = cell(order, 1);

            % Initialize indices for FoW
            obj.indicesFoW = cell(order, 1);
            for o = 2:order
                obj.indicesFoW{o} = obj.get_indices_FoW(o);
            end
        end
        
        function obj = compute_manifold(obj)
            % Compute the LSM manifold.

            % Order 1
            obj.R{1} = sparse(obj.sys.Lambda);
            obj.gamma{1} = obj.sys.omega0; % 1/2 * imag([1, -1] * sum(obj.R{1}, 2));
            obj.W{1} = sparse(obj.sys.V);

            % Higher orders
            for ii = 2:obj.order
                % Find C matrix
                obj.C{ii} = obj.compute_C(ii);

                % Normal form parametrization
                obj.R{ii} = obj.normal_form_parametrization(ii);

                % Gamma values: gamma = 1/2 * imag([1, -1] * sum(obj.R{ii}, 2));
                obj.gamma{ii} = full(imag(sum(obj.R{ii}(1, :))));

                % Cohomological equation
                obj.W{ii} = obj.solve_cohomological(ii);
            end
        end

        % Compute Y matrix
        function Y = compute_Y(obj, ii, jj)
            % Compute Y matrix for the computation of C matrix.
            % Inputs:
            %   ii: row index (current order).
            %   jj: column index (previous orders).
            % Outputs:
            %   Y: the Y matrix.

            Y = sparse(obj.Q^jj, obj.Q^ii);
            thisR = obj.R{ii - jj + 1};
            if nnz(thisR) > 0
                for kk = 1:jj
                    % Left identity
                    if kk > 1
                        leftI = speye(obj.Q^(kk - 1));
                    else
                        leftI = 1;
                    end
                    
                    % Right identity
                    if kk < jj
                        rightI = speye(obj.Q^(jj - kk));
                    else
                        rightI = 1;
                    end

                    % Compute Y
                    Y = Y + kron(kron(leftI, thisR), rightI); % kronn({leftI, thisR, rightI});
                end
            end
        end

        % Compute C matrix
        function C = compute_C(obj, order)
            % Compute the C matrix for the computation of the W matrix.
            % Inputs:
            %   order: the order of the C matrix.
            % Outputs:
            %   C: the C matrix.

            C = obj.compute_FoW(order);
            for jj = 2:(order - 1)
                C = C - obj.sys.B * obj.W{jj} * obj.compute_Y(order, jj);
            end
        end

        function indices = get_indices_FoW(obj, order)
            % Get indices for the computation of (F o W) based on the order.
            % Inputs:
            %   order: the current order.
            % Outputs:
            %   indices: the indices for the computation of (F o W) at the current order.

            indices = cell(0);
            for ii = 2:length(obj.sys.F) % terms to sum
                v = (1:(order - ii + 1)).';
                data = disp_with_rep(v, ii);
                indices = [indices, data(sum(data, 2) == order, :)];
            end
        end

        function FoW = compute_FoW(obj, order)
            % Compute (F o W) matrix.
            % Inputs:
            %   order: the order of the (F o W) matrix.
            % Outputs:
            %   FoW: the (F o W) matrix.

            % Extract indices
            indices = obj.indicesFoW{order};

            % Loop 
            FoW = sparse(obj.sys.N, obj.Q^order);
            for ii = 1:length(indices)
                idx = indices{ii};
                cols = size(idx, 2);
                if nnz(obj.sys.F{cols}) > 0
                    for jj = 1:size(idx, 1)
                        FoW = FoW + obj.sys.F{cols} * kronn(obj.W(idx(jj, :)));
                    end
                end
            end
        end

        function [E, N] = compute_EN(obj, order)
            % Compute E and N matrices for the computation of the R matrix.
            % These matrices depends on the inner resonances.
            % Inputs:
            %   order: the order of the E and N matrices.
            % Outputs:
            %   E: the E matrix.
            %   N: the N matrix.

            % Find inner resonances
            [lIndices, jIndices] = find_inner_resonance(order);
            nInnerResonances = length(lIndices);

            % Initialize E and N triplets
            tripletsE = zeros(nInnerResonances, 3); % each column of E has nnz = 1
            tripletsN = zeros(nInnerResonances * obj.sys.N, 3); % each column of N has nnz = obj.sys.N

            % Loop over inner resonances
            itE = 0; itN = 0;
            for ii = 1:nInnerResonances
                % Triplets of E
                tripletsE(itE + 1, :) = [obj.Q*(lIndices(ii) - 1) + jIndices(ii), ii, 1];
                itE = itE + 1;

                % Triplets of N
                tripletsN(itN + (1:obj.sys.N), :) = [obj.sys.N*(lIndices(ii) - 1) + (1:obj.sys.N).', ii*ones(obj.sys.N, 1), obj.sys.U(:, jIndices(ii))];
                itN = itN + obj.sys.N;
            end

            % Build sparse matrices from triplets
            E = sparse(tripletsE(:, 1), tripletsE(:, 2), tripletsE(:, 3), obj.Q^order * obj.Q, nInnerResonances);
            N = sparse(tripletsN(:, 1), tripletsN(:, 2), tripletsN(:, 3), obj.Q^order * obj.sys.N, nInnerResonances);
        end

        function R = normal_form_parametrization(obj, order)
            % Compute the normal form parametrization coefficients.
            % Inputs:
            %   order: the order of the R matrix.
            % Outputs:
            %   R: the R matrix that contains the normal form parametrization coefficients.

            % Check order
            if mod(order, 2) == 0 % for even orders R = 0
                R = sparse(obj.Q, obj.Q^order);
            else
                [E, N] = obj.compute_EN(order);
                Rvec = E * N' * obj.C{order}(:);
                R = reshape(Rvec, obj.Q, []);
            end
        end

        function [L, Y] = compute_L(obj, order)
            % Compute L matrix.
            % Inputs:
            %   order: the order of the L matrix.
            % Outputs:
            %   L: the L matrix.
            %   Y: the Y matrix.

            Y = obj.compute_Y(order, order);
            L = kron(Y.', obj.sys.B) - kron(speye(obj.Q^order), obj.sys.A);
        end

        function D = compute_D(obj, order)
            % Compute D matrix
            % Inputs:
            %   order: the order of the D matrix.
            % Outputs:
            %   D: the D matrix.

            D = kron(speye(obj.Q^order), obj.sys.B*obj.W{1});
        end

        function [h, D] = compute_h(obj, order)
            % Compute h vector.
            % Inputs:
            %   order: the order of the h vector.
            % Outputs:
            %   h: the h vector.
            %   D: the D matrix.

            D = obj.compute_D(order);
            h = obj.C{order}(:) - D*obj.R{order}(:);
        end

        function W = solve_cohomological(obj, order)
            % Solve the cohomological equation.
            % Inputs:
            %   order: the order of the W matrix.
            % Outputs:
            %   W: the W matrix.

            % Find L and h
            L = obj.compute_L(order);
            h = obj.compute_h(order);

            % Solve
            Wvec = lsqminnorm(L, h);
            W = reshape(Wvec, obj.sys.N, []);
        end

        function omega = compute_omega(obj, rho)
            % Compute omega given a rho value.
            % Inputs:
            %   rho: the rho values.
            % Outputs:
            %   omega: the omega values.

            omega = zeros(size(rho));
            for ii = 1:obj.order
                omega = omega + obj.gamma{ii} .* rho.^(ii - 1);
            end
        end
        
        function [z, zk] = compute_z_single(obj, rho, iDof)
            % Compute the rms value of z for a single rho value.
            % Inputs:
            %   rho: the rho value.
            %   iDof: the degree of freedom index.
            % Outputs:
            %   z: the rms value of z.
            %   zk: the z values for each theta.

            % Check input
            if numel(rho) ~= 1
                error('The input rho should have only one element.');
            end

            % Loop over theta
            zk = zeros(1, obj.nTheta);
            z = 0;
            for kk = 1:obj.nTheta

                % Define the reduced coordinates
                pk = rho * exp(obj.theta(kk).*[1j; -1j]);

                % First order
                pkI = pk;
                zk(kk) = real(obj.W{1}(iDof, :) * pkI);

                % Loop over orders to build zk
                for ii = 2:obj.order
                    pkI = kron(pkI, pk);
                    zk(kk) = zk(kk) + real(obj.W{ii}(iDof, :) * pkI);
                end

                % Update z
                z = z + zk(kk).^2;
            end

            % Compute RMS
            z = sqrt(z ./ obj.nTheta);
        end

        function z = compute_z(obj, rho, iDof)
            % Compute the rms value of z for a range of rho values.
            % Inputs:
            %   rho: the rho values.
            %   iDof: the degree of freedom index.
            % Outputs:
            %   z: the rms values of z.

            % Loop over rho
            z = zeros(size(rho));
            for jj = 1:length(rho)
                z(jj) = obj.compute_z_single(rho(jj), iDof);
            end
        end

        function [c, ceq] = nonlcon(obj, rho, iDof, z)
            % Nonlinear equality constraints for fmincon.
            c = [];
            ceq = obj.compute_z_single(rho, iDof) - z;
        end

        function rho = solve_rho_at_leading_order(obj, z, iDof)
            % Solve for the rho value at the leading order.
            % Inputs:
            %   z: the rms values of z.
            %   iDof: the degree of freedom index.
            % Outputs:
            %   rho: the rho values.

            % Loop over z
            rho = zeros(size(z));
            w1 = full(obj.W{1}(iDof, :));
            for jj = 1:length(z)
                zkTildeQuadSum = 0;
                for kk = 1:obj.nTheta
                    pk_tilde = exp(obj.theta(kk).*[1j; -1j]);
                    zkTildeQuadSum = zkTildeQuadSum + real(w1 * pk_tilde).^2;
                end
                rho(jj) = z(jj) ./ sqrt(zkTildeQuadSum / obj.nTheta);
            end            
        end

        % function rho = solve_rho(obj, z, iDof)
        %     % Solve for the rho values given the rms values of z.
        %     % Inputs:
        %     %   z: the rms values of z.
        %     %   iDof: the degree of freedom index.
        %     % Outputs:
        %     %   rho: the rho values.
        % 
        %     % Rho step
        %     drho = obj.solve_rho_at_leading_order(max(z) / 1000, iDof);
        % 
        %     % Build the (rho, z) pairs
        %     ii = 1;
        %     rhoTemp = 0;
        %     zTemp = 0;
        %     r = 0;
        %     while r < 1
        %         rhoTemp(ii + 1) = rhoTemp(ii) + drho;
        %         zTemp(ii + 1) = obj.compute_z_single(rhoTemp(ii + 1), iDof);
        %         ii = ii + 1;
        %         r = max(zTemp) / max(z);
        %     end
        % 
        %     % Interpolation
        %     rho = zeros(size(z));
        %     for ii = 1 : length(z)
        %         [~, ind] = min(abs(zTemp - z(ii)));
        %         rho(ii) = rhoTemp(ind);
        %     end
        % end

        function rho = solve_rho(obj, z, iDof)
            % Solve for the rho values given the rms values of z.
            % Inputs:
            %   z: the rms values of z.
            %   iDof: the degree of freedom index.
            % Outputs:
            %   rho: the rho values.

            % Initialize rho
            rho = obj.solve_rho_at_leading_order(z, iDof);

            % Loop over z
            for ii = 1:length(z)
                if z(ii) == 0
                    rho(ii) = 0;
                else
                    % Solve with fzero
                    fun = @(rho) obj.compute_z_single(rho, iDof) - z(ii);
                    rho(ii) = fzero(fun, rho(ii));
                end
            end            
        end
    end
end
