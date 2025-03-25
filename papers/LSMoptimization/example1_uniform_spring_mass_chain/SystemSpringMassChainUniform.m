classdef SystemSpringMassChainUniform < SystemBase
    % A class that describes the uniform spring-mass chain system.
    % Methods:
    %   SystemSpringMassChainUniform: constructor.
    %   second_order_form: initialize the second order form.
    %   compute_nonlinear_triplets: compute the nonlinear force triplets.

    properties
        nMasses, nSprings, bc
        dofs_to_dofs_reduced
        dofs_reduced_to_dofs
        M_tilde, K_tilde, f_tilde
    end
    
    methods
        % Constructor
        function obj = SystemSpringMassChainUniform(x, args)
            % Constructor.
            % Inputs:
            %   x: design variables.
            %   args: additional arguments.

            % Extract variables
            m = x(1);
            k = x(2);
            k2 = x(3);
            k3 = x(4);

            % Number of design variables
            obj.ndv = length(x);

            % Extract auxiliar variables
            obj.bc = args.bc;

            % Second order form
            obj = obj.second_order_form(m, k, k2, k3);

            % State space form z = [x, dx]
            obj = obj.state_space();
        end

        % Initialize second order form
        function obj = second_order_form(obj, m, k, k2, k3)
            % Initialize the second order form.
            % Inputs:
            %   m: mass.
            %   k: linear stiffness.
            %   k2: quadratic stiffness.
            %   k3: cubic stiffness.
            
            % Number of masses and springs
            obj.nMasses = length(obj.bc); % total number of masses (free and fixed)
            obj.nSprings = obj.nMasses - 1;
            obj.n = sum(obj.bc == 1); % number of free masses

            % Map between reduced and fixed dofs
            obj.dofs_to_dofs_reduced = -ones(obj.nMasses, 1);
            obj.dofs_to_dofs_reduced(obj.bc == 1) = 1:obj.n;
            temp = 1:obj.nMasses;
            obj.dofs_reduced_to_dofs = temp(obj.bc == 1);

            % Build mass matrix
            M_tr = zeros(obj.n, 3);
            jj = 0;
            for ii = 1:obj.nMasses
                if obj.dofs_to_dofs_reduced(ii) > -1
                    jj = jj + 1;
                    M_tr(jj, :) = [obj.dofs_to_dofs_reduced(ii), obj.dofs_to_dofs_reduced(ii), 1];
                end
            end
            obj.M_tilde = sparse(M_tr(:, 1), M_tr(:, 2), M_tr(:, 3), obj.n, obj.n);
            obj.M = sparse(M_tr(:, 1), M_tr(:, 2), M_tr(:, 3) * m, obj.n, obj.n);

            % Build stiffness matrix
            K_tr = -ones(4*obj.nSprings, 3);
            Ke_tr = [0, 0,  1;
                     1, 0, -1;
                     0, 1, -1;
                     1, 1,  1];
            kk = 0;
            for ii = 1:obj.nSprings
                for jj = 1:size(Ke_tr, 1)
                    row = ii + Ke_tr(jj, 1);
                    col = ii + Ke_tr(jj, 2);
                    val = Ke_tr(jj, 3);
                    if obj.dofs_to_dofs_reduced(row) > -1 && obj.dofs_to_dofs_reduced(col) > -1
                        kk = kk + 1;
                        K_tr(kk, :) = [obj.dofs_to_dofs_reduced(row), obj.dofs_to_dofs_reduced(col), val];
                    end
                end
            end
            obj.K_tilde = sparse(K_tr(1:kk, 1), K_tr(1:kk, 2), K_tr(1:kk, 3), obj.n, obj.n);
            obj.K = sparse(K_tr(1:kk, 1), K_tr(1:kk, 2), K_tr(1:kk, 3) * k, obj.n, obj.n);

            % Nonlinear force vectors
            obj.f_tilde = cell(3, 1);
            obj.f = cell(3, 1);

            % Order 1
            obj.f_tilde{1} = sparse(obj.n, 2*obj.n);
            obj.f{1} = sparse(obj.n, 2*obj.n);

            % Higher orders
            k_nl = [0, k2, k3];
            for ii = 2:length(k_nl)
                [rows, cols, vals] = obj.compute_nonlinear_triplets(ii);
                obj.f_tilde{ii} = sparse(rows, cols, vals, obj.n, (2*obj.n)^ii);
                obj.f{ii} = sparse(rows, cols, vals * k_nl(ii), obj.n, (2*obj.n)^ii);
            end
        end

        function [rows, cols, vals] = compute_nonlinear_triplets(obj, order)
            % Compute the nonlinear force triplets for a given order.
            % Inputs:
            %   order: order of the nonlinear force.
            % Outputs:
            %   rows: row indices of the triplets.
            %   cols: column indices of the triplets.
            %   vals: values of the triplets.
            
            % Initialize triplets
            rows = []; cols = []; vals = [];

            % Initialize indices
            idx = disp_with_rep([0; 1], order);
        
            % Loop over dofs
            for dof = 1:(obj.nMasses - 1)
        
                % Loop over the coefficients of the polynomial
                for jj = 1:size(idx, 1)
            
                    % Extract indices of the x values involved in this coefficient
                    idx_loc = idx(jj, :);
                    idx_glob = dof + idx_loc;
            
                    % Check if all the indices are not fixed
                    if sum(obj.bc(idx_glob)) == length(idx_glob)
            
                        % Compute the global index of the coefficient
                        idx_red = obj.dofs_to_dofs_reduced(idx_glob);
                        col = (2*obj.n).^(order - (1:order)) * (idx_red(:) - 1) + 1;
            
                        % Add to the triplet at the row of the left dof
                        if obj.dofs_to_dofs_reduced(dof) > -1
                            rows = [rows; obj.dofs_to_dofs_reduced(dof)];
                            cols = [cols; col];
                            vals = [vals; (-1)^(sum(idx_loc) + order - 1)]; % this is the value of the nonlinear force
                        end
            
                        % Add to the triplet at the row of the right dof
                        if obj.dofs_to_dofs_reduced(dof + 1) > -1
                            rows = [rows; obj.dofs_to_dofs_reduced(dof + 1)];
                            cols = [cols; col];
                            vals = [vals; (-1)^(sum(idx_loc) + order)]; % this is the value of the nonlinear force
                        end
                    end
                end
            end
        end

        function obj = second_order_form_sensitivity(obj)
            % Initialize the sensitivity of the second order form.

            % M matrix
            obj.dM = cell(1, obj.ndv);
            obj.dM{1} = obj.M_tilde;
            obj.dM{2} = sparse(obj.n, obj.n);
            obj.dM{3} = sparse(obj.n, obj.n);
            obj.dM{4} = sparse(obj.n, obj.n);

            % K matrix
            obj.dK = cell(1, obj.ndv);
            obj.dK{1} = sparse(obj.n, obj.n);
            obj.dK{2} = obj.K_tilde;
            obj.dK{3} = sparse(obj.n, obj.n);
            obj.dK{4} = sparse(obj.n, obj.n);

            % f matrix
            obj.df = cell(length(obj.f), obj.ndv);

            % f1 matrix
            obj.df{1, 1} = sparse(obj.n, 2*obj.n);
            obj.df{1, 2} = sparse(obj.n, 2*obj.n);
            obj.df{1, 3} = sparse(obj.n, 2*obj.n);
            obj.df{1, 4} = sparse(obj.n, 2*obj.n);

            % f2 matrix
            obj.df{2, 1} = sparse(obj.n, (2*obj.n)^2);
            obj.df{2, 2} = sparse(obj.n, (2*obj.n)^2);
            obj.df{2, 3} = obj.f_tilde{2};
            obj.df{2, 4} = sparse(obj.n, (2*obj.n)^2);

            % f3 matrix
            obj.df{3, 1} = sparse(obj.n, (2*obj.n)^3);
            obj.df{3, 2} = sparse(obj.n, (2*obj.n)^3);
            obj.df{3, 3} = sparse(obj.n, (2*obj.n)^3);
            obj.df{3, 4} = obj.f_tilde{3};
        end
    end
end
