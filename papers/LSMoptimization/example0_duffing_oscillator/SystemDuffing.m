classdef SystemDuffing < SystemBase
    % A class that describes the Duffing oscillator system.
    % Methods:
    %   SystemDuffing: constructor.
    %   second_order_form: initialize the second order form.
    %   second_order_form_sensitivity: initialize the sensitivity.
    
    methods
        function obj = SystemDuffing(x, args)
            % Constructor.
            % Inputs:
            %   x: design variables.
            %   args: additional arguments.

            % Check inputs
            if nargin > 1
                args = [];
            end

            % Extract variables
            m = x(1);
            k = x(2);
            k2 = x(3);
            k3 = x(4);
            obj.ndv = 4;

            % Second order form
            obj = obj.second_order_form(m, k, k2, k3);

            % State space form z = [x, dx]
            obj = obj.state_space();
        end

        function obj = second_order_form(obj, m, k, k2, k3)
            % Initialize the second order form.
            % Inputs:
            %   m: mass.
            %   k: stiffness.
            %   k2: stiffness.
            %   k3: stiffness.

            obj.n = 1;
            obj.M = sparse(1, 1, m, obj.n, obj.n);
            obj.K = sparse(1, 1, k, obj.n, obj.n);
            obj.f{1} = sparse(obj.n, 2*obj.n);
            obj.f{2} = sparse(1, 1, k2, obj.n, (2*obj.n)^2);
            obj.f{3} = sparse(1, 1, k3, obj.n, (2*obj.n)^3);
        end

        function obj = second_order_form_sensitivity(obj)
            % Initialize the sensitivity of the second order form.

            % M matrix
            obj.dM = cell(1, obj.ndv);
            obj.dM{1} = sparse(1, 1, 1, obj.n, obj.n);
            obj.dM{2} = sparse(obj.n, obj.n);
            obj.dM{3} = sparse(obj.n, obj.n);
            obj.dM{4} = sparse(obj.n, obj.n);

            % K matrix
            obj.dK = cell(1, obj.ndv);
            obj.dK{1} = sparse(obj.n, obj.n);
            obj.dK{2} = sparse(1, 1, 1, obj.n, obj.n);
            obj.dK{3} = sparse(obj.n, obj.n);
            obj.dK{4} = sparse(obj.n, obj.n);

            % f1 matrix
            obj.df{1, 1} = sparse(obj.n, 2*obj.n);
            obj.df{1, 2} = sparse(obj.n, 2*obj.n);
            obj.df{1, 3} = sparse(obj.n, 2*obj.n);
            obj.df{1, 4} = sparse(obj.n, 2*obj.n);

            % f2 matrix
            obj.df{2, 1} = sparse(obj.n, (2*obj.n)^2);
            obj.df{2, 2} = sparse(obj.n, (2*obj.n)^2);
            obj.df{2, 3} = sparse(1, 1, 1, obj.n, (2*obj.n)^2);
            obj.df{2, 4} = sparse(obj.n, (2*obj.n)^2);

            % f3 matrix
            obj.df{3, 1} = sparse(obj.n, (2*obj.n)^3);
            obj.df{3, 2} = sparse(obj.n, (2*obj.n)^3);
            obj.df{3, 3} = sparse(obj.n, (2*obj.n)^3);
            obj.df{3, 4} = sparse(1, 1, 1, obj.n, (2*obj.n)^3);
        end
    end
end
