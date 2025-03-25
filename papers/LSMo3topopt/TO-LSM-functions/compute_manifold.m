function [Omega, theta, wManifold, wBB] = compute_manifold(M, K, T2, T3, omega, phi, rho, idx, doPlot)
    % Compute the 3rd order manifold and the relative backbone curve
    % Inputs:
    %   M: mass matrix
    %   K: stiffness matrix
    %   T2: 2nd order tensor
    %   T3: 3rd order tensor
    %   omega: frequency
    %   phi: mode shape
    %   rhoBB: amplitude in the reduced space
    %   idx: index of the degrees of freedom
    %   doPlot: flag to plot the results (optional, default is false)
    % Outputs:
    %   Omega: frequency
    %   theta: angles
    %   wManifold: manifold in the physical space
    %   wBB: backbone curve in the physical space

    % Check inputs
    if nargin < 9
        doPlot = false;
    end

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

    % Backbone in the reduced space
    Omega = omega + gamma * rho.^2;

    % Theta points
    nTheta = 2^7;
    theta = linspace(0, 2*pi, nTheta);

    % Compute manifold
    wManifold = zeros(length(rho), length(theta));
    for ii = 1:length(rho)
        for kk = 1:length(theta)
            p = rho(ii) * exp([1i; -1i] * theta(kk));
            wManifold(ii, kk) = phi(idx) * p(1) + phi(idx) * p(2) + ...
                w20(idx) * p(1)^2 + w11(idx) * p(1) * p(2) + w02(idx) * p(2)^2 + ...
                w30(idx) * p(1)^3 + w21(idx) * p(1)^2 * p(2) + ...
                w12(idx) * p(1) * p(2)^2 + w03(idx) * p(2)^3;
        end
    end
    wManifold = real(wManifold);

    % Compute backbone in the physical space
    wBB = sqrt(sum(wManifold.^2, 2) / nTheta);

    % Plot the results
    if doPlot
        % Plot manifold
        figure
        hold on; grid on; box on;
        surf(rho.' * cos(theta), rho.' * sin(theta), wManifold, ...
            'FaceColor', [0.0745 0.6235 1], 'EdgeColor', 'none', 'FaceAlpha', 0.7)
        L = light;
        L.Position = [0 0 1];
        lighting gouraud
        rhoIndex = round(length(rho)/3 * [1, 2]);
        for ii = rhoIndex
            x = rho(ii) * cos(theta); y = rho(ii) * sin(theta);
            plot3(x, y, wManifold(ii, :), 'r.-', 'LineWidth', 2)
        end
        xlabel('$\rho \cos(\theta)$', 'Interpreter', 'latex')
        ylabel('$\rho \sin(\theta)$', 'Interpreter', 'latex')
        title('Manifold', 'Interpreter', 'latex')
        zlabel('$w$', 'Interpreter', 'latex')
        set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
        axis square tight
        box on
        view(45,15)

        % Plot backbone
        figure
        hold on; grid on;
        plot(Omega, wBB, 'LineWidth', 2)
        xlabel('$\Omega$', 'Interpreter', 'latex')
        ylabel('$w$', 'Interpreter', 'latex')
        title('Backbone curve', 'Interpreter', 'latex')
        set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
        axis square tight
        box on
    end
end
