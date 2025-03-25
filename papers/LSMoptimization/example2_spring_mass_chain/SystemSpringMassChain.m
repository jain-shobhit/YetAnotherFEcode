classdef SystemSpringMassChain < SystemBase
    % A class that describes the uniform spring-mass chain system.
    % Methods:
    %   SystemSpringMassChain: constructor.
    %   second_order_form: initialize the second order form.
    %   compute_nonlinear_triplets_single_spring: compute the nonlinear force triplets for a single spring.
    %   compute_nonlinear_triplets: compute the nonlinear force triplets.
    %   second_order_form_sensitivity: initialize the sensitivity of the second order form.

    properties
        nMasses, nSprings, bc
        dofs_to_dofs_reduced
        dofs_reduced_to_dofs
    end
    
    methods
        function obj = SystemSpringMassChain(x, args)
            % Constructor.
            % Inputs:
            %   x: design variables.
            %   args: additional arguments.

            % Size of the chain
            obj.bc = args.bc;
            obj.nMasses = length(args.bc);
            obj.nSprings = obj.nMasses - 1;

            % Number of design variables
            obj.ndv = length(x);

            % Second order form
            obj = obj.second_order_form(x);

            % State space form z = [x, dx]
            obj = obj.state_space();
        end

        function obj = second_order_form(obj, x)
            % Initialize the second order form.
            % Inputs:
            %   x: design variables.

            % Number of dof
            obj.n = sum(obj.bc == 1);

            % Extract properties
            m = zeros(size(obj.bc)); % all the masses (free and fixed)
            m(obj.bc == 1) = x(1:obj.n);
            k = x(obj.n + (1:obj.nSprings));
            k2 = x(obj.n + obj.nSprings + (1:obj.nSprings));
            k3 = x(obj.n + 2*obj.nSprings + (1:obj.nSprings));

            % Map between reduced and fixed dofs
            obj.dofs_to_dofs_reduced = -ones(obj.nMasses, 1);
            obj.dofs_to_dofs_reduced(obj.bc == 1) = 1:obj.n;
            temp = 1:obj.nMasses;
            obj.dofs_reduced_to_dofs = temp(obj.bc == 1);

            % Build mass matrix
            tripletsM = zeros(obj.n, 3);
            jj = 0;
            for ii = 1:obj.nMasses
                if obj.dofs_to_dofs_reduced(ii) > -1
                    jj = jj + 1;
                    tripletsM(jj, :) = [obj.dofs_to_dofs_reduced(ii), obj.dofs_to_dofs_reduced(ii), m(ii)];
                end
            end
            obj.M = sparse(tripletsM(:, 1), tripletsM(:, 2), tripletsM(:, 3), obj.n, obj.n);

            % Build stiffness matrix
            tripletsK = zeros(4*obj.nSprings, 3);
            tripletsKe = [0, 0,  1;
                          1, 0, -1;
                          0, 1, -1;
                          1, 1,  1]; % element matrix of one spring
            kk = 0;
            for ii = 1:obj.nSprings % loop over springs
                for jj = 1:size(tripletsKe, 1)
                    row = ii + tripletsKe(jj, 1);
                    col = ii + tripletsKe(jj, 2);
                    val = tripletsKe(jj, 3) * k(ii);
                    if obj.dofs_to_dofs_reduced(row) > -1 && obj.dofs_to_dofs_reduced(col) > -1
                        kk = kk + 1;
                        tripletsK(kk, :) = [obj.dofs_to_dofs_reduced(row), obj.dofs_to_dofs_reduced(col), val];
                    end
                end
            end
            obj.K = sparse(tripletsK(1:kk, 1), tripletsK(1:kk, 2), tripletsK(1:kk, 3), obj.n, obj.n);

            % Nonlinear force vectors
            obj.f = cell(3, 1);

            % Order 1
            obj.f{1} = sparse(obj.n, 2*obj.n);

            % Higher orders
            kNl = [zeros(obj.nSprings, 1), k2(:), k3(:)];
            for ii = 2:size(kNl, 2)
                [rows, cols, vals] = obj.compute_nonlinear_triplets(ii, kNl(:, ii));
                obj.f{ii} = sparse(rows, cols, vals, obj.n, (2*obj.n)^ii);
            end
        end

        function [rows, cols, vals] = compute_nonlinear_triplets_single_spring(obj, order, idxK, valK)
            % Compute the nonlinear force triplets for a single spring.
            % Inputs:
            %   order: order of the polynomial spring.
            %   idxK: index of the spring.
            %   valK: value of the spring.
            % Outputs:
            %   rows: row indices of the triplets.
            %   cols: column indices of the triplets.
            %   vals: values of the triplets.

            % Initialize variables
            idx = disp_with_rep([0; 1], order);
            rows = []; cols = []; vals = [];
        
            % Loop over the coefficients of the polynomial
            for jj = 1:size(idx, 1)
                % Extract indices of the x values involved in this coefficient
                idxLocal = idx(jj, :);
                idxGlobal = idxK + idxLocal;
        
                % Check if all the indices are not fixed
                if sum(obj.bc(idxGlobal)) == length(idxGlobal)
                    % Compute the global index of the coefficient
                    idx_red = obj.dofs_to_dofs_reduced(idxGlobal);
                    col = (2*obj.n).^(order - (1:order)) * (idx_red(:) - 1) + 1;
        
                    % Add to the triplet at the row of the left mass
                    if obj.dofs_to_dofs_reduced(idxK) > -1
                        rows = [rows; obj.dofs_to_dofs_reduced(idxK)];
                        cols = [cols; col];
                        vals = [vals; (-1)^(sum(idxLocal) + order - 1) * valK]; % this is the value of the nonlinear force
                    end
        
                    % Add to the triplet at the row of the right mass
                    if obj.dofs_to_dofs_reduced(idxK + 1) > -1
                        rows = [rows; obj.dofs_to_dofs_reduced(idxK + 1)];
                        cols = [cols; col];
                        vals = [vals; (-1)^(sum(idxLocal) + order) * valK]; % this is the value of the nonlinear force
                    end
                end
            end
        end

        function [rows, cols, vals] = compute_nonlinear_triplets(obj, order, k)
            % Compute the nonlinear force triplets for a given order.
            % Inputs:
            %   order: order of the nonlinear force.
            %   k: stiffness of the springs.
            % Outputs:
            %   rows: row indices of the triplets.
            %   cols: column indices of the triplets.
            %   vals: values of the triplets.
            
            % Initialize triplets
            rows = []; cols = []; vals = [];
        
            % Loop over springs
            for ii = 1:obj.nSprings
                % Compute the triplets associated to the current spring
                [row, col, val] = compute_nonlinear_triplets_single_spring(obj, order, ii, k(ii));

                % Update triplets
                rows = [rows; row];
                cols = [cols; col];
                vals = [vals; val];
            end
        end

        function obj = second_order_form_sensitivity(obj)
            % Initialize the sensitivity of the second order form.
            
            % Design variables
            obj.ndv = obj.n + 3 * obj.nSprings;

            % M matrix
            obj.dM = cell(1, obj.ndv);
            for ii = 1:obj.ndv
                if ii <= obj.n
                    obj.dM{ii} = sparse(ii, ii, 1, obj.n, obj.n);
                else
                    obj.dM{ii} = sparse(obj.n, obj.n);
                end
            end

            % K matrix
            obj.dK = cell(1, obj.ndv);
            for ii = 1:obj.n
                obj.dK{ii} = sparse(obj.n, obj.n); % derivative wrt the masses
            end
            tripletsK = zeros(4, 3);
            tripletsKe = [0, 0,  1;
                          1, 0, -1;
                          0, 1, -1;
                          1, 1,  1]; % element matrix of one spring (triplet form)
            for ii = 1:obj.nSprings % loop over springs

                % Reset counter
                kk = 0;

                % Loop over the elements of tripletsKe
                for jj = 1:size(tripletsKe, 1)

                    % Global indices and derivative value
                    row = ii + tripletsKe(jj, 1); % global row
                    col = ii + tripletsKe(jj, 2); % global col

                    % Check if the current mass is free
                    row_reduced = obj.dofs_to_dofs_reduced(row);
                    col_reduced = obj.dofs_to_dofs_reduced(col);
                    if row_reduced > -1 && col_reduced > -1

                        % Update counter and store triplet
                        kk = kk + 1; 
                        tripletsK(kk, :) = [row_reduced, col_reduced, tripletsKe(jj, 3)];
                    end
                end

                % Derivative wrt the linear springs
                obj.dK{obj.n + ii} = sparse(tripletsK(1:kk, 1), tripletsK(1:kk, 2), tripletsK(1:kk, 3), obj.n, obj.n);
            end
            for ii = (obj.n + obj.nSprings + 1):obj.ndv
                obj.dK{ii} = sparse(obj.n, obj.n); % derivative wrt the nonlinear springs
            end

            % Derivative of the nonlinear terms
            for ii = 1:obj.ndv
                obj.df{1, ii} = sparse(obj.n, 2*obj.n);
            end
            springIndex = 1:obj.nSprings;
            for ii = 2:length(obj.f) % ii represents the order

                % The derivatives wrt to masses and linear springs are null
                for jj = 1:(obj.n + obj.nSprings)
                    obj.df{ii, jj} = sparse(obj.n, (2*obj.n)^ii);
                end

                % Find the global indices of the current order springs
                springIndexGlobal = obj.n + (ii - 1)*obj.nSprings + springIndex;

                % Loop over nonlinear springs
                for jj = (obj.n + obj.nSprings + 1):obj.ndv % jj is the global index of design variables

                    % Check if the jj-th spring belongs to the current order
                    if sum(springIndexGlobal == jj) == 1

                        % Local spring index
                        springIndexLocal = springIndex(springIndexGlobal == jj);
    
                        % Compute the triplets for the current spring at this order
                        [rows, cols, vals] = compute_nonlinear_triplets_single_spring(obj, ii, springIndexLocal, 1);

                        % Assemble sparse matrix
                        obj.df{ii, jj} = sparse(rows, cols, vals, obj.n, (2*obj.n)^ii);
                    else

                        % Null sensitivity if the jj-th spring does not belongs to the current order 
                        obj.df{ii, jj} = sparse(obj.n, (2*obj.n)^ii);
                    end
                end
            end
        end
    end
end