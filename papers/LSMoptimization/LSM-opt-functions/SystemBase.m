classdef SystemBase
    % A base class for a mechanical system used in the optimization.
    % Properties:
    %   n: the number of degrees of freedom.
    %   M: the mass matrix.
    %   K: the stiffness matrix.
    %   f: cell array with the nonlinear force terms.
    %   omega: vector of eigenfrequencies.
    %   phi: matrix of mode shapes (n x nev).
    %   omega0: the target eigenfrequency.
    %   phi0: the target mode shape.
    %   N: the size of the state vector.
    %   A: the state matrix.
    %   B: the input matrix.
    %   F: cell array with the nonlinear force terms.
    %   Lambda: the state eigenvalue matrix (2 x 2).
    %   V: the state right eigenvector matrix (N x nev).
    %   U: the state left eigenvector matrix (N x nev).
    %   ndv: the number of design variables.
    %   dM: cell array with the mass matrix sensitivities.
    %   dK: cell array with the stiffness matrix sensitivities.
    %   df: cell matrix with the nonlinear force term sensitivities.
    %   dA: cell array with the state matrix sensitivities.
    %   dB: cell array with the input matrix sensitivities.
    %   dF: cell matrix with the nonlinear force term sensitivities.
    % Methods:
    %   modal_analysis: compute the eigenfrequencies and mode shapes of the second order system.
    %   state_space: compute the state space form of the system.
    %   state_space_sensitivity: compute the sensitivity of the state space form of the system.
    %   modal_analysis_sensitivity: compute the sensitivity of the target eigenfrequency and mode shape.
    %   state_space_modal_analysis: compute the eigenvalues and eigenvectors of the state space form of the system.
    %   state_space_modal_analysis_sensitivity: compute the sensitivity of the eigenvalues and eigenvectors of the state space form of the system.
    %   compute_nonlinear_force: compute the nonlinear force vector.
    %   compute_numerical_frequency_response: compute the frequency response of the system using numerical integration.
    %   frequency_sweep: performs negative and positive frequency sweeps.

    properties
        % Second order system
        n
        M, K, f
        omega, phi
        omega0 = NaN, phi0 = NaN

        % State space form
        N
        A, B, F
        Lambda, V, U
        
        % Sensitivity
        ndv
        dM, dK, df
        dA, dB, dF
    end

    methods
        function obj = modal_analysis(obj, nev)
            % Compute the eigenfrequencies and mode shapes (mass-normalized) of the system.
            % Inputs:
            %   nev: the number of eigenvalues to compute. Default is 10 or n, whichever is smaller.

            % Number of modes
            if nargin < 2
                nev = min(10, obj.n);
            end

            % Solve the inverted eigenvalue problem
            [Phi, D] = eigs(obj.M, obj.K, nev, 'largestabs');

            % Eigenfrequencies
            obj.omega = sqrt(1 ./ diag(D));

            % Mass normalization
            obj.phi = zeros(size(Phi));
            for ii = 1:size(Phi, 2)
                xMx = Phi(:, ii).' * obj.M * Phi(:, ii);
                obj.phi(:, ii) = Phi(:, ii) / sqrt(xMx);
            end
        end

        function obj = state_space(obj)
            % Compute the state space form of the system.
            % The state space form is given by:
            %   B = [-K, 0; 0, M]
            %   A = [0, -K; -K, 0]
            %   F = [0; -f]
            %
            obj.N = 2 * obj.n;
            obj.B = [-obj.K, sparse(obj.n, obj.n); sparse(obj.n, obj.n), obj.M];
            obj.A = [sparse(obj.n, obj.n), -obj.K; -obj.K, sparse(obj.n, obj.n)];
            for ii = 1:length(obj.f)
                obj.F{ii} = [sparse(obj.n, obj.N^ii); -obj.f{ii}];
            end
        end

        function obj = state_space_sensitivity(obj)
            % Compute the sensitivity of the state space form of the system.

            % Initialize
            obj.dB = cell(1, obj.ndv);
            obj.dA = cell(1, obj.ndv);
            obj.dF = cell(length(obj.F), obj.ndv);

            % Loop over design variables
            for dv = 1:obj.ndv
                % Matrices B and A
                obj.dB{dv} = [-obj.dK{dv}, sparse(obj.n, obj.n); sparse(obj.n, obj.n), obj.dM{dv}];
                obj.dA{dv} = [sparse(obj.n, obj.n), -obj.dK{dv}; -obj.dK{dv}, sparse(obj.n, obj.n)];

                % Nonlinear terms
                for ii = 1:length(obj.f)
                    obj.dF{ii, dv} = [sparse(obj.n, obj.N^ii); -obj.df{ii, dv}];
                end
            end
        end

        function [dphi0, domega0] = modal_analysis_sensitivity(obj)
            % Compute the sensitivity of the target eigenfrequency and mode shape.
            % Outputs:
            %   dphi0: the sensitivity of the target mode shape.
            %   domega0: the sensitivity of the target eigenfrequency.

            % Check if target eigenfrequency is defined
            if isnan(obj.omega0)
                error('Target eigenfrequency is not defined.')
            end

            % Initialize sensitivity
            domega0 = cell(1, obj.ndv);
            dphi0 = cell(1, obj.ndv);

            % Build matrix
            mat = [obj.K - obj.omega0^2*obj.M, -2*obj.omega0*obj.M*obj.phi0;
                   -2*obj.omega0*obj.phi0.'*obj.M, 0];

            % Build right-hand sides
            rhs = zeros(obj.n + 1, obj.ndv);
            for ii = 1:obj.ndv
                rhs(:, ii) = [-(obj.dK{ii} - obj.omega0^2*obj.dM{ii})*obj.phi0;
                              obj.omega0*obj.phi0.'*obj.dM{ii}*obj.phi0];
            end

            % Solve linear system
            X = mat\rhs;

            % Extract sensitivity
            for ii = 1:obj.ndv
                dphi0{ii} = X(1:end-1, ii);
                domega0{ii} = X(end, ii);
            end
        end

        function obj = state_space_modal_analysis(obj)
            % Computes the igenvalues and eigenvectors of the state space form of the system.
            % Lambda = omega0 * diag([1j, -1j])
            % V = [phi0, phi0; 1j*omega0*phi0, -1j*omega0*phi0]
            % U = -1/(2*omega0^2) * conj(V)

            obj.Lambda = obj.omega0 * diag([1j, -1j]);
            obj.V = [obj.phi0, obj.phi0; 1j*obj.omega0*obj.phi0, -1j*obj.omega0*obj.phi0];
            obj.U = -1/(2*obj.omega0^2) * conj(obj.V);
        end

        function [LambdaTilde, VTilde, dVTilde, UTilde, dUTilde] = state_space_modal_analysis_sensitivity(obj)
            % Computes the sensitivity of the eigenvalues and eigenvectors of the state space form of the system.
            % Outputs:
            %   LambdaTilde: normalized version of the state eigenvalue matrix.
            %   VTilde: normalized version of the state right eigenvector matrix.
            %   dVTilde: sensitivity of VTilde with respect to the eigenfrequency.
            %   UTilde: normalized version of the state left eigenvector matrix.
            %   dUTilde: sensitivity of UTilde with respect to the eigenfrequency.
            
            LambdaTilde = sparse(diag([1j, -1j]));
            VTilde = sparse([1, 1; 1j*obj.omega0, -1j*obj.omega0]);
            dVTilde = sparse([0, 0; 1j, -1j]);
            UTilde = -1/(2*obj.omega0^2) * conj(VTilde);
            dUTilde = sparse([1/obj.omega0^3, 1/obj.omega0^3; -1j/(2*obj.omega0^2), 1j/(2*obj.omega0^2)]);
        end

        function fNl = compute_nonlinear_force(obj, z)
            % Compute the nonlinear force vector.
            % Inputs:
            %   z: the state vector.
            % Outputs:
            %   Fnl: the nonlinear force vector.

            zk = z;
            fNl = zeros(obj.N, 1);
            for ii = 2:length(obj.F)
                zk = kron(zk, z);
                if nnz(obj.F{ii}) > 0
                    fNl = fNl + obj.F{ii} * zk;
                end
            end
        end

        function [omegaNeg, omegaPos, zNeg, zPos] = compute_numerical_frequency_response(obj, omegaLim, nOmega, fExtNorm, iInput, iOutput, rayleighDamping, doPrint)
            % Compute the frequency response of the system using numerical integration.
            % Inputs:
            %   omegaLim: the frequency limits (ratio of the eigenfrequency).
            %   nOmega: the number of frequency points.
            %   fExtNorm: vector of magnitude of the external force.
            %   iInput: the index of the input degree of freedom.
            %   iOutput: the index of the output degree of freedom.
            %   rayleighDamping: the Rayleigh damping coefficients (alpha, beta).
            %   doPrint: the print flag. Default is false.
            % Outputs:
            %   omegaNeg: the negative sweep frequency vector.
            %   omegaPos: the positive sweep frequency vector.
            %   zNeg: the negative sweep response.
            %   zPos: the positive sweep response.

            % Damping matrix
            C = rayleighDamping(1) * obj.M + rayleighDamping(2) * obj.K;
            stateC = [sparse(obj.n, obj.n), sparse(obj.n, obj.n);
                      sparse(obj.n, obj.n), -C];

            % Q-factor and transient time
            csi = rayleighDamping(1) / (2*obj.omega0) + rayleighDamping(2) * obj.omega0 / 2;
            if csi~=0
                Q = 1/(2*csi);
            else
                Q = 0; % to avoid inf transient
            end
            
            % Exponential decay from linear theory:
            %     u = u0*exp(-pi*f0*t/Q) 
            %     log(u/u0) = -pi*f0*t/Q
            %     t = -Q/(pi*f0)*log(u/u0)
            % finally, setting u/u0 = 5%:
            f0 = obj.omega0 / 2 / pi;
            tTransient = -Q/(pi*f0)*log(0.05);

            % State matrix
            dampedA = obj.A + stateC;

            % Positive and negative frequency sweeps
            omegaLim = omegaLim * obj.omega0;
            omegaNeg = linspace(omegaLim(2), omegaLim(1), nOmega);
            omegaPos = linspace(omegaLim(1), omegaLim(2), nOmega);

            % External force
            fExt = zeros(obj.N, 1);

            % Ode options
            options = odeset('RelTol', 1e-4, 'AbsTol', 1e-4, 'Mass', obj.B);

            % Loop over amplitude values
            zNeg = cell(length(fExtNorm), 1);
            zPos = cell(length(fExtNorm), 1);
            for ii = 1:length(fExtNorm)

                % Display message
                if doPrint
                    tic;
                    fprintf('Forcing amplitude %d ... ', ii)
                end
                
                % External force
                fExt(iInput) = fExtNorm(ii);

                % Negative sweep
                if doPrint
                    fprintf('negative sweep ... ')
                end
                zNeg{ii} = obj.frequency_sweep(dampedA, omegaNeg, fExt, iOutput, tTransient, options);

                % Positive sweep
                if doPrint
                    fprintf('positive sweep ... ')
                end
                zPos{ii} = obj.frequency_sweep(dampedA, omegaPos, fExt, iOutput, tTransient, options);

                % Final message
                if doPrint
                    dt = toc;
                    fprintf('done in %.4f seconds.\n', dt)
                end
            end
        end

        function z = frequency_sweep(obj, dampedA, omega, fExt, iDof, tTransient, options)
            % Performs negative and positive frequency sweeps.
            % Inputs:
            %   dampedA: the damped state matrix.
            %   omega: the frequency vector.
            %   fExt: the external force vector.
            %   iDof: the degree of freedom index.
            %   tTransient: the transient time.
            %   options: the ode options.
            % Outputs:
            %   z: the response vector.
            
            % Loop over omega
            nCycles = 50; % number of cycles
            nPpC = 25; % number of points per cycle
            z = zeros(size(omega)); X0 = zeros(obj.N, 1);
            for ii = 1:length(omega)

                % Time vector
                T = 2*pi / omega(ii);
                nTransientCycles = ceil(tTransient / T);
                t = linspace(0, (nCycles + nTransientCycles)*T, nPpC * (nCycles + nTransientCycles));

                % Equation of motion
                odefun = @(t, X) dampedA * X + obj.compute_nonlinear_force(X) + fExt * cos(omega(ii) * t);

                % Solve ode
                [~, X] = ode23(odefun, t, X0, options);
                
                % Update initial conditions with the last final conditions
                X0 = X(end, :).';

                % Extract the trajectory of the target mass and compute rms
                x = X(t > tTransient, iDof);
                z(ii) = rms(x - mean(x)); % the mean is used to remove constant components
            end
        end
    end
end