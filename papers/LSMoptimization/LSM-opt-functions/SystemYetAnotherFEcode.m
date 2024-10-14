classdef SystemYetAnotherFEcode < SystemBase
    % A class that describes a von Karman beam defined using Yet Another FE Code (YAFEC).
    % Properties:
    %   yafec_assembly_constructor: function that produces the yafec assembly having as input a parameter vector p.
    %   myAssembly: YAFEC model.
    %   nodes: nodes of the YAFEC model.
    %   T2, T3: stiffness tensors associated to quadratic and cubic terms.
    %   F2, F3: stiffness matrices associated to quadratic and cubic terms (1st order form).
    % Methods:
    %   SystemYetAnotherFEcode: constructor.
    %   second_order_form: initialize the second order form.
    %   second_order_form_sensitivity: initialize the sensitivity of the second order form.
    %   tensor2matrix_1stOrd: converts the stiffness tensors into their corresponding matrix representation (1st order form).

    properties
        yafec_assembly_constructor
        myAssembly
        nodes
        T2, T3
        F2, F3
        p
    end

    methods
        function obj = SystemYetAnotherFEcode(p, args)
            % Constructor.
            % Inputs:
            %   p: design variables.
            %   args: additional arguments.
            %   args.yafec_assembly_constructor: function that produces the
            %   yafec assembly having as input a parameter vector p.

            % Check inputs
            if nargin == 1
                args = [];
            end

            % Define the YAFEC model
            obj.yafec_assembly_constructor = args.yafec_assembly_constructor;
            obj.myAssembly = feval(obj.yafec_assembly_constructor, p);

            % Nodes
            obj.nodes = obj.myAssembly.Mesh.nodes;

            % Current parameter vector
            obj.p = p;
            obj.ndv = length(obj.p);

            % Second order form
            obj = obj.second_order_form();

            % State space form z = [x, dx]
            obj = obj.state_space();
        end

        function obj = second_order_form(obj)
            % Initialize the second order form.

            % Extract the number of degrees of freedom
            n = obj.myAssembly.Mesh.nDOFs;

            % Initial zero displacement
            u0 = zeros(n, 1);

            % Mass and stiffness matrices
            M = obj.myAssembly.mass_matrix();
            K = obj.myAssembly.tangent_stiffness_and_force(u0);

            % Stiffness tensors
            t2 = obj.myAssembly.tensor('T2', [n n n], [2, 3]);
            t3 = obj.myAssembly.tensor('T3', [n n n n], [2, 3, 4]);

            % Constrained matrices and tensors
            obj.K = obj.myAssembly.constrain_matrix(K);
            obj.M = obj.myAssembly.constrain_matrix(M);
            obj.T2 = obj.myAssembly.constrain_tensor(t2);
            obj.T3 = obj.myAssembly.constrain_tensor(t3);

            % Number of constrained dofs
            obj.n = size(obj.K, 1);

            % Convert stiffness tensors to matrix representation
            % obj = obj.tensor2matrix_2ndOrd();
            obj = obj.tensor2matrix_1stOrd();

            % Remove zeros (they are re-added by default in the SystemBase.state_space() method)
            obj.f{1} = sparse(obj.n, 2 * obj.n);
            obj.f{2} = obj.F2(1:obj.n, :);
            obj.f{3} = obj.F3(1:obj.n, :);
            % note: in SystemBase.StateSpaceForm(obj), F will come with a minus sign
        end

        function obj = second_order_form_sensitivity(obj)
            % Initialize the sensitivity of the second order form.
            
            % Perturbation and nominal parameter vector
            dp = 1e-4;
            p0 = obj.p;

            % Loop over design variables
            args.yafec_assembly_constructor = obj.yafec_assembly_constructor;
            for ii = 1 : obj.ndv
                % Perturbed parameter vector
                pNew = p0;
                pNew(ii) = pNew(ii) + 1e-4;

                % Perturbed object
                objNew = SystemYetAnotherFEcode(pNew, args);

                % Finite differences (Forward Euler)
                obj.dM{ii} = (objNew.M - obj.M) / dp;
                obj.dK{ii} = (objNew.K - obj.K) / dp;
                obj.df{1, ii} = sparse(obj.n, 2*obj.n);
                obj.df{2, ii} = (objNew.f{2} - obj.f{2}) / dp;
                obj.df{3, ii} = (objNew.f{3} - obj.f{3}) / dp;
            end
        end

        function obj = tensor2matrix_1stOrd(obj)
            % Converts the stiffness tensors into their corresponding matrix representation.
            % For example, being z the state vector:
            %   quadratic nonlinear force: ttv(T2, {u, u}, [2, 3])       <--> F2 * kronexp(z, 2)
            %   cubic nonlinear force:     ttv(T3, {u, u, u}, [2, 3, 4]) <--> F3 * kronexp(z, 3)

            % x: displacement vector, size(x) = [n, 1]
            % z: state vector, z = [x; \dot{x}], size(z) = [2n, 1]
            N = 2 * obj.n;

            % Concatenate tensors
            T23 = {obj.T2, obj.T3};

            % Loop over tensors
            for o = 3:4
                T = T23{o - 2};
                % Loop over 1st dimension of T
                nEl = nnz(T); % nnz elements of tensor T
                triplets = zeros(nEl, 3); % triplet list [rows, cols, vals]
                for ii = 1:size(T.subs, 1)
                    % Row index
                    row = T.subs(ii, 1);

                    % Column index of the short matrix f(x)
                    col = 1;
                    for jj = o:-1:2
                        % col = col + obj.n^(o - jj) * (T.subs(ii, jj) - 1); % 2nd order matrix
                        col = col + N^(o - jj) * (T.subs(ii, jj) - 1); % 1st order matrix
                    end

                    % Value
                    triplets(ii, :) = [row, col, T.vals(ii)];
                end
                if o == 3
                    obj.F2 = sparse(triplets(:, 1), triplets(:, 2), triplets(:, 3), N, N^(o - 1));
                else
                    obj.F3 = sparse(triplets(:, 1), triplets(:, 2), triplets(:, 3), N, N^(o - 1));
                end
            end
        end
    end
end