function plot_manifold(lsm, rhoTarget, iDof)
    % Plot the manifold of the Spectral Submanifold object.
    % Inputs:
    %   lsm: the Lyapunov Subcenter Manifolds object.
    %   rhoTarget: vector, the target amplitudes in the reduced space.
    %   iDof: scalar, the degree of freedom index.

    % Number of points
    nTheta = 61;
    nRho = 15;

    % Meshgrid
    rho = linspace(0, 1.2 * max(rhoTarget), nRho);
    theta = linspace(0, 2*pi, nTheta);

    % Compute physical coordinates for all rho and theta
    z = zeros(nRho, nTheta);
    for ii = 1 : nRho
        for jj = 1 : nTheta
            p = rho(ii) * [exp( 1i*theta(jj));
                        exp(-1i*theta(jj))];
            z(ii, jj) = eval_z(lsm, p, iDof);
        end
    end

    % Plot manifold
    hold on; grid on; box on;
    surf(rho.' .* cos(theta), rho.' .* sin(theta), z, ...
        'FaceColor',[0.0745 0.6235 1], 'EdgeColor','none','FaceAlpha',0.7)
    L = light;
    L.Position = [0 0 1];
    lighting gouraud
    axis square tight
    xlabel('$\rho \cos(\theta)$', 'Interpreter', 'latex')
    ylabel('$\rho \sin(\theta)$', 'Interpreter', 'latex')
    zlabel('$Z$ [m]', 'Interpreter', 'latex')
    set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
    view(45,15)

    % Plot orbits
    for ii = 1:length(rhoTarget)
        % Compute orbit
        z = zeros(1, nTheta);
        for jj = 1 : nTheta
            p = rhoTarget(ii) * exp([ 1i*theta(jj);
                                     -1i*theta(jj)]);
            z(jj) = eval_z(lsm, p, iDof);
        end

        % Plot
        plot3(rhoTarget(ii) * cos(theta), rhoTarget(ii) * sin(theta), z, 'r-', 'LineWidth', 1.5)
    end
end

% Auxiliary function used to evaluate the manifold
function z = eval_z(lsm, p, iDof)
    % Evaluate physical coordinates of the manifold given the reduced coordinates.
    % Inputs:
    %   lsm: the Lyapunov Subcenter Manifolds object.
    %   p: the reduced coordinates.
    %   iDof: the degree of freedom index.
    % Outputs:
    %   z: the physical coordinates.

    % First order
    pExp = p;
    z = lsm.W{1}(iDof, :) * pExp;

    % Loop over higher orders
    for ii = 2:lsm.order
        pExp = kron(pExp, p);
        z = z + lsm.W{ii}(iDof, :) * pExp;
    end

    % Real part
    z = real(z);
end
