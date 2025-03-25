function rho = compute_rho(w, M, K, T2, T3, omega, phi, idx, rho0)
    % Compute the reduced displacement (rho) given the physical one (w)
    % Inputs:
    %   w: physical displacement
    %   M: mass matrix
    %   K: stiffness matrix
    %   T2: 2nd order tensor
    %   T3: 3rd order tensor
    %   omega: frequency
    %   phi: mode shape
    %   idx: index of the degrees of freedom
    %   rho0: initial guess for the reduced displacement (optional)
    % Outputs:
    %   rho: reduced displacement

    % Order 2: index 20
    Lambda20 = 2i * omega;
    f20 = ttv(T2, {phi, phi}, [3, 2]);
    f20 = double(f20);
    L20 = K + Lambda20^2 * M;
    w20 = lsqminnorm(L20, -f20);
    
    % Order 2: index 11
    f11 = 2 * f20;
    L11 = K;
    w11 = lsqminnorm(L11, -f11);
    
    % Order 2: index 02
    w02 = conj(w20);

    % Order 3: index 30
    Lambda30 = 3i * omega;
    f30 = ttv(T3, {phi, phi, phi}, [4, 3, 2]) + ttv(T2, {w20, phi}, [3, 2]) + ttv(T2, {w20, phi}, [2, 3]);
    f30 = double(f30);
    L30 = K + Lambda30^2 * M;
    w30 = lsqminnorm(L30, -f30);

    % Order 3: index 21
    Lambda21 = 1i * omega;
    f21 = ttv(T2, {w20 + w11, phi}, [3, 2]) + ...
          ttv(T2, {w20 + w11, phi}, [2, 3]) + ...
          3 * ttv(T3, {phi, phi, phi}, [4, 3, 2]);
    f21 = double(f21);
    gamma = 1/(2*omega) * phi.' * f21;
    h21 = f21 - M * (2*omega) * phi * gamma;
    L21 = K + Lambda21^2 * M;
    w21 = lsqminnorm(L21, -h21);
    
    % Order 2: index 02
    w12 = conj(w21);
    
    % Order 2: index 03
    w03 = conj(w30);

    % Theta points
    nTheta = 2^7;
    theta = linspace(0, 2*pi, nTheta);

    % Solve with fzero
    if w == 0
        rho = 0;
    else
        fun = @(x) compute_w(x, phi, w20, w11, w02, w30, w21, w12, w03, idx, theta) - w;
        if nargin < 9
            rho0 = compute_rho_leading_order(w, phi, idx, theta);
        end
        rho = fzero(fun, rho0);
    end
end

% Helper functions
function w = compute_w(rho, phi, w20, w11, w02, w30, w21, w12, w03, idx, theta)
    % Compute the physical displacement given the reduced one
    % Inputs:
    %   rho: reduced displacement
    %   phi: mode shape
    %   w20: manifold coefficient for index 20
    %   w11: manifold coefficient for index 11
    %   w02: manifold coefficient for index 02
    %   w30: manifold coefficient for index 30
    %   w21: manifold coefficient for index 21
    %   w12: manifold coefficient for index 12
    %   w03: manifold coefficient for index 03
    %   idx: index of the degrees of freedom
    %   theta: vector of angles
    % Outputs:
    %   w: physical displacement
    
    % Compute manifold
    w = 0;
    nTheta = length(theta);
    for kk = 1:nTheta
        p = rho * exp([1i; -1i] * theta(kk));
        wManifold = phi(idx) * p(1) + phi(idx) * p(2) + ...
            w20(idx) * p(1)^2 + w11(idx) * p(1) * p(2) + w02(idx) * p(2)^2 + ...
            w30(idx) * p(1)^3 + w21(idx) * p(1)^2 * p(2) + ...
            w12(idx) * p(1) * p(2)^2 + w03(idx) * p(2)^3;
        w = w + real(wManifold).^2;
    end
    w = sqrt(w / nTheta);
end

function rho = compute_rho_leading_order(w, phi, idx, theta)
    % Compute the reduced displacement (rho) given the physical one (w)
    % at the leading order (only index 10 and 01)
    % Inputs:
    %   w: physical displacement
    %   phi: mode shape
    %   idx: index of the degrees of freedom
    %   theta: vector of angles
    % Outputs:
    %   rho: reduced displacement
    
    % Theta points
    nTheta = length(theta);

    % Compute manifold
    wTilde = 0;
    for kk = 1:nTheta
        p = exp([1i; -1i] * theta(kk));
        wManifold = phi(idx) * p(1) + phi(idx) * p(2);
        wTilde = wTilde + real(wManifold).^2;
    end
    rho = w / sqrt(wTilde / nTheta);
end
