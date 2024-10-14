function [x, orderNew] = run_fmincon(x0, A, b, Aeq, beq, xl, xu, sysArgs, lsmArgs, optArgs, options)
    % Run the optimization with 'fmincon'.
    % This allows to use global variables to update the expansion order.
    % Inputs:
    %   x0: the initial design variables.
    %   A: the linear inequality constraint matrix.
    %   b: the linear inequality constraint vector.
    %   Aeq: the linear equality constraint matrix.
    %   beq: the linear equality constraint vector.
    %   xl: the lower bound vector.
    %   xu: the upper bound vector.
    %   sysArgs: the system arguments structure.
    %   lsmArgs: the LSM arguments structure.
    %   optArgs: the optimization arguments structure.
    %   options: the optimization options.
    % Outputs:
    %   x: the optimal design variables.
    %   orderNew: the new expansion order.

    % Define output function
    options.OutputFcn = @(x, optimValues, state) output_fun(x, optimValues, state, sysArgs, optArgs);

    % Solve with fmincon
    x = fmincon(@(x) fun(x), x0, A, b, Aeq, beq, xl, xu, @(x) nonlcon(x, sysArgs, optArgs), options);

    % Update expansion order
    orderNew = lsmArgs.order;

    %% Objective function and sensitivities
    function [fval, fgrad] = fun(x)
        fval = x(2) * x(4);
        fgrad = [0, x(4), 0, x(2)];
    end
    
    %% Nonlinear constraint function and sensitivities
    function [c, ceq, dc, dceq] = nonlcon(x, sysArgs, optArgs)
        % Nonlinear constraints function for optimization. This function computes
        % the equality and inequality constraints and their sensitivities.
        % Inputs:
        %   x: the design variables.
        %   sysArgs: the arguments structure.
        %   optArgs: the optimization arguments structure.
        % Outputs:
        %   c: the inequality constraints.
        %   ceq: the equality constraints.
        %   dc: the inequality constraints sensitivities.
        %   dceq: the equality constraints sensitivities.
        
        % Build system
        sys = feval(sysArgs.systemType, x, sysArgs);
        
        % Modal analysis
        sys = sys.modal_analysis();
        
        % Identify master mode
        if isfield(sysArgs, 'modeIndex')
            % Use the provided mode index
            modeIndex = sysArgs.modeIndex;
        elseif isfield(sysArgs, 'phiRef')
            % Use the modal assurance criterion
            modeIndex = modal_assurance_criterion(sys.phi, sysArgs.phiRef);
        else
            % Use the first mode (default)
            modeIndex = 1;
        end
        
        % Assign master mode
        sys.omega0 = sys.omega(modeIndex);
        sys.phi0 = sys.phi(:, modeIndex);
        sys = sys.state_space_modal_analysis();
        
        % Compute LSM
        lsm = LSM(sys, lsmArgs.order);
        lsm = lsm.compute_manifold();
        
        % Compute LSM sensitivities
        sens = Sensitivity(sys, lsm);
        sens = sens.sensitivity_lsm();
        
        % Exctract nonlinear constraints
        iDof = lsmArgs.iDof;
        zTarget = optArgs.zTarget;
        omegaTarget = optArgs.omegaTarget;
        coeffTarget = optArgs.coeffTarget;
        
        % Find number of equality and inequality constraints
        iCeq = find(coeffTarget == 0);
        iC = find(abs(coeffTarget) == 1);
        
        % Resize vectors
        ceq = zeros(1, length(iCeq));
        dceq = zeros(length(x), length(iCeq));
        c = zeros(1, length(iC));
        dc = zeros(length(x), length(iC));
        
        % Loop over equality constraints
        for ii = 1:length(iCeq)
            [fval, fgrad] = sens.sensitivity_backbone(zTarget(iCeq(ii)), iDof);
            ceq(ii) = fval - omegaTarget(iCeq(ii));
            dceq(:, ii) = fgrad;
        end
        
        % Loop over inequality constraints
        for ii = 1:length(iC)
            [fval, fgrad] = sens.sensitivity_backbone(zTarget(iC(ii)), iDof);
            c(ii) = (fval - omegaTarget(iC(ii))) * coeffTarget(iC(ii));
            dc(:, ii) = fgrad * coeffTarget(iC(ii));
        end
    end
    
    %% Output function for fmincon
    function stop = output_fun(x, optimValues, state, sysArgs, optArgs)
        % Output function for the optimization with 'fmincon'.
        % Inputs:
        %   x: the design variables.
        %   optimValues: the optimization values.
        %   state: the optimization state (init, iter, done).
        %   sysArgs: the system arguments structure.
        %   optArgs: the optimization arguments structure.
        % Outputs:
        %   stop: the stop flag.
        
        % Output
        stop = false;
        
        % Iteration number
        iter = optimValues.iteration;
        
        % Build system
        sys = feval(sysArgs.systemType, x, sysArgs);
        
        % Modal analysis
        sys = sys.modal_analysis();
        
        % Identify master mode
        if isfield(sysArgs, 'modeIndex')
            % Use the provided mode index
            modeIndex = sysArgs.modeIndex;
        elseif isfield(sysArgs, 'phiRef')
            % Use the modal assurance criterion
            modeIndex = modal_assurance_criterion(sys.phi, sysArgs.phiRef);
        else
            % Use the first mode (default)
            modeIndex = 1;
        end
        
        % Assign master mode
        sys.omega0 = sys.omega(modeIndex);
        sys.phi0 = sys.phi(:, modeIndex);
        sys = sys.state_space_modal_analysis();
        
        % Compute LSM
        lsm = LSM(sys, lsmArgs.order);
        lsm = lsm.compute_manifold();
        
        % Compute omega
        rhoTarget = lsm.solve_rho(optArgs.zTarget, lsmArgs.iDof);
        rho = linspace(0, max(rhoTarget), 51);
        z = lsm.compute_z(rho, lsmArgs.iDof);
        omega = lsm.compute_omega(rho);
        
        % Switch state
        switch state
            case 'init'
                % Initialize backbone plot
                figure();
                hold on; grid on; box on;
                plot_backbone(omega, z, optArgs.omegaTarget, optArgs.zTarget)
                axis tight;
                xlabel('$\Omega$ [rad/s]', 'Interpreter', 'latex')
                ylabel('$Z$ [m]', 'Interpreter', 'latex')
                drawnow;
            case 'iter'
                % Evaluate the residual error of the invariance equation
                err = error_measure(lsm, lsmArgs, optArgs);
        
                % Tolerance on the residual error (adjust this value if needed)
                tol = 0.018;
        
                % Check the residual error
                warning('off', 'backtrace')
                if err > tol && lsmArgs.order < lsmArgs.orderMax
                    lsmArgs.order = lsmArgs.order + 2;
                    warning('Residual is %.3f. Order increased (%d->%d). Continuing.', err, lsmArgs.order - 2, lsmArgs.order)
                elseif err > tol && lsmArgs.order >= lsmArgs.orderMax
                    warning('Residual is %.3f but max order has been reached. Continuing.', err)
                elseif err < tol / 100
                    lsmArgs.order = lsmArgs.order - 1;
                    warning('Residual is %.3f. Order reduced (%d->%d). Continuing.', err, lsmArgs.order + 1, lsmArgs.order)    
                end
        
                % Display iteration and design variables
                disp(['Iter ', num2str(iter), ': ', num2str(x), ' - res: ', num2str(err)]);
        
                % Plot backbone
                plot_backbone(omega, z, optArgs.omegaTarget, optArgs.zTarget)
                drawnow;
            case 'done'
                % Display final iteration and design variables
                disp(['Done: ', num2str(x)]);
        
                % Plot backbone
                plot_backbone(omega, z, optArgs.omegaTarget, optArgs.zTarget)
                drawnow;
            otherwise
        end
    end
end