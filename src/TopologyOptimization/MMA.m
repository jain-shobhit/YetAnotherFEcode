classdef MMA < handle
    % A class to handle the Method of Moving Asymptotes (MMA).
    % The original code is available under the
    % GNU General Public License (GPLv3) at: https://www.smoptit.se/
    % The source files are includesd in "external/GCMMA-MMA-code-1.5".
    % Properties:
    %   m, n: the number of constraints and design variables.
    %   xmin, xmax: the lower and upper bounds for the variables x_j.
    %   xold1, xold2: xval, one and two iterations ago.
    %   low, upp: the lower and upper asymptotes from the previous
    %       iteration.
    %   a0, a, c, d: the constants in the terms a_0*z, a_i*z, c_i*y_i,
    %       0.5*d_i*(y_i)^2.
    %   move: maximum change in the design variables.
    %   freeElements: free elements indicator: 1 free, 0 fixed (optional,
    %       default is []).
    %   freeElementsFlag: flag for the free elements indicator (optional,
    %       default is []).
    %   mapSymmetry: the SymmetryMap object (optional, default is []).
    %   mapSymmetryFlag: flag for symmetry map (optional, default is []).
    % Methods:
    %   MMA: constructor.
    %   optimize: perform the MMA optimization.
    
    properties
        m, n
        xmin, xmax
        xold1, xold2
        low, upp
        a0, a, c, d
        move
        freeElements = []
        freeElementsFlag = false
        mapSymmetry = []
        mapSymmetryFlag = false
    end
    
    methods
        function obj = MMA(m, move, to, varargin)
            % MMA: constructor.
            % Inputs:
            %   m: the number of general constraints.
            %   move: maximum change in the design variables.
            %   to: the TopologyOptimization object.
            %   mapSymmetry: SymmetryMap object (optional, default is []).
            % Outputs:
            %   obj: the MMA object.

            % Parse inputs
            p = inputParser;
            addOptional(p, 'mapSymmetry', []);
            parse(p, varargin{:});
            obj.mapSymmetry = p.Results.mapSymmetry;

            % Store inputs
            obj.freeElements = to.freeElements;
            x = to.d;

            % Check symmetry
            if ~isempty(obj.mapSymmetry)
                obj.freeElements = obj.mapSymmetry.symmetry_extraction(obj.freeElements);
                x = obj.mapSymmetry.symmetry_extraction(x);
                obj.mapSymmetryFlag = true;
            end

            % Check free elements
            obj.freeElements = (obj.freeElements == 1);
            if ~all(obj.freeElements)
                % obj.freeElements = (obj.freeElements == 1);
                x = x(obj.freeElements);
                obj.freeElementsFlag = true;
            end

            % Number of design variables
            n = length(x);
            
            % Initialize
            obj.m     = m;                  % The number of general constraints.
            obj.n     = n;                  % The number of design variables x_j.
            obj.xmin  = zeros(n,1);         % Column vector with the lower bounds for the variables x_j.
            obj.xmax  = ones(n,1);          % Column vector with the upper bounds for the variables x_j.
            obj.xold1 = x(:);               % xval, one iteration ago (provided that iter>1).
            obj.xold2 = x(:);               % xval, two iterations ago (provided that iter>2).
            obj.low   = ones(n,1);          % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
            obj.upp   = ones(n,1);          % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).
            obj.a0    = 1;                  % The constants a_0 in the term a_0*z.
            obj.a     = zeros(m,1);         % Column vector with the constants a_i in the terms a_i*z.
            obj.c     = 10000*ones(m,1);    % Column vector with the constants c_i in the terms c_i*y_i.
            obj.d     = zeros(m,1);         % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
            obj.move = move;                % Maximum change in the design variables.
        end

        function x = optimize(obj, iter, xval, f0val, df0dx, fval, dfdx)
            % Perform the MMA optimization.
            % Inputs:
            %   iter: the number of iterations.
            %   xval: the design variables.
            %   f0val: the objective function value.
            %   df0dx: the gradient of the objective function.
            %   fval: the constraint function values.
            %   dfdx: the gradients of the constraint functions.
            % Outputs:
            %   x: the optimized design variables.

            % Apply symmetry reduction
            if obj.mapSymmetryFlag
                xval = obj.mapSymmetry.symmetry_extraction(xval);
                df0dx = obj.mapSymmetry.symmetry_reduction(df0dx);
                dfdx = obj.mapSymmetry.symmetry_reduction(dfdx.').';
            end
            
            % Remove fixed elements
            if obj.freeElementsFlag
                x0 = xval;
                xval = xval(obj.freeElements);
                df0dx = df0dx(obj.freeElements);
                dfdx = dfdx(:, obj.freeElements);
            end

            % MMA step
            [x, ~, ~, ~, ~, ~, ~, ~, ~, low_new, upp_new] = ...
                mmasub(obj.m, obj.n, iter, xval, obj.xmin, obj.xmax, obj.xold1, obj.xold2, ...
                       f0val,df0dx,fval,dfdx,obj.low,obj.upp,obj.a0,obj.a,obj.c,obj.d, obj.move);

            % MMA update
            obj.xold2 = obj.xold1;
            obj.xold1 = xval;
            obj.low = low_new;
            obj.upp = upp_new;

            % Update design variables (only free elements)
            if obj.freeElementsFlag
                x0(obj.freeElements) = x;
                x = x0;
            end

            % Apply symmetry extension
            if obj.mapSymmetryFlag
                x = obj.mapSymmetry.symmetry_extension(x);
            end
        end
    end
end
