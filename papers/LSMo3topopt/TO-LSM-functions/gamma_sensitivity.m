function [gamma, dgamma] = gamma_sensitivity(myAssembly, M, K, T2, T3, omega, phi, Me, Ke, T2e, T3e, sensPhK, sensPhM, freeElements)
    % Compute the 3rd order backbone coefficient (gamma) and its adjoint sensitivity (dgamma) for a topology optimization problem
    % Inputs:
    %   myAssembly: instance of the Assembly class
    %   M: mass matrix
    %   K: stiffness matrix
    %   T2: 2nd order tensor
    %   T3: 3rd order tensor
    %   omega: frequency
    %   phi: mode shape
    %   Me: element mass matrix
    %   Ke: element stiffness matrix
    %   T2e: element 2nd order tensor
    %   T3e: element 3rd order tensor
    %   sensPhK: physical density sensitivity of the stiffness matrix
    %   sensPhM: physical density sensitivity of the mass matrix
    %   freeElements: flag to indicate which elements are free
    % Outputs:
    %   gamma: 3rd order backbone coefficient
    %   dgamma: adjoint sensitivity of the 3rd order backbone coefficient

    % Check if the physical density sensitivities are provided
    if nargin == 11
        sensPhK = ones(myAssembly.Mesh.nElements, 1);
        sensPhM = ones(myAssembly.Mesh.nElements, 1);
        freeElements = ones(myAssembly.Mesh.nElements, 1);
    elseif nargin == 13
        freeElements = ones(myAssembly.Mesh.nElements, 1);
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
    
    % Compute gamma
    f21 = ttv(T2, {w20 + w11, phi}, [3, 2]) + ...
          ttv(T2, {w20 + w11, phi}, [2, 3]) + ...
          3 * ttv(T3, {phi, phi, phi}, [4, 3, 2]);
    f21 = double(f21);
    gamma = 1/(2*omega) * phi.' * f21;

    % Partial derivative wrt w20
    f21_w20 = ttv(T2, {phi}, 2) + ttv(T2, {phi}, 3);
    gamma_w20 = 1/(2*omega) * double(ttv(f21_w20, {phi}, 1));
    
    % Partial derivative wrt w11
    gamma_w11 = gamma_w20;
    
    % Partial derivative wrt omega
    gamma_omega = -gamma / omega;
    
    % Partial derivative wrt phi ( = phi = phi)
    f20_phi = f21_w20;
    w20_w11 = w20 + w11;
    f21_phi = ttv(T2, {w20_w11}, 2) + ttv(T2, {w20_w11}, 3) + ...
              3 * ttv(T3, {phi, phi}, [3, 2]) + ...
              3 * ttv(T3, {phi, phi}, [4, 2]) + ...
              3 * ttv(T3, {phi, phi}, [4, 3]);
    gamma_phi = 1/(2*omega) * (f21 + double(ttv(f21_phi, {phi}, 1)));

    % Adjoint of w20
    lambda20 = lsqminnorm(L20.', -gamma_w20);
    
    % Adjoint of w11
    lambda11 = lsqminnorm(L11.', -gamma_w11);
    
    % Adjoint of phi and omega
    b_omega = gamma_omega - 8 * omega * lambda20.' * M * w20;
    b_phi = gamma_phi + double(ttv(f20_phi, lambda20 + 2 * lambda11, 1));
    b = [-b_phi; b_omega / omega];
    A = [(K - omega^2 * M), 2 * M * phi;
         2 * phi.' * M, 0];
    x = A\b;
    lambda0 = x(1:end-1);
    lambda1 = x(end);

    % Unconstrained vectors
    phiUnc = myAssembly.unconstrain_vector(phi);
    w20Unc = myAssembly.unconstrain_vector(w20);
    w11Unc = myAssembly.unconstrain_vector(w11);
    lambda0Unc = myAssembly.unconstrain_vector(lambda0);
    lambda20Unc = myAssembly.unconstrain_vector(lambda20);
    lambda11Unc = myAssembly.unconstrain_vector(lambda11);

    % Index of free elements
    freeElementsIdx = find(freeElements == 1).';

    % Sensitivity
    dgamma = zeros(myAssembly.Mesh.nElements, 1);
    for ii = freeElementsIdx
        % Dofs of the current element
        iDOFs = myAssembly.Mesh.Elements(ii).Object.iDOFs;

        % Extract element quantities
        phiE = phiUnc(iDOFs, 1);
        w20E = w20Unc(iDOFs, 1);
        w11E = w11Unc(iDOFs, 1);
        lambda0E = lambda0Unc(iDOFs, 1);
        lambda20E = lambda20Unc(iDOFs, 1);
        lambda11E = lambda11Unc(iDOFs, 1);

        % Physical density sensitivities
        dK = sensPhK(ii) * Ke;
        dM = sensPhM(ii) * Me;
        dT2 = sensPhK(ii) * T2e;
        dT3 = sensPhK(ii) * T3e;

        % Matrices and tensors derivatives
        dL20 = dK + Lambda20^2 * dM;
        dL11 = dK;
        dT2dT2 = dT2 + permute(dT2, [1, 3, 2]);

        % Partial derivatives
        f20E = double(ttv(dT2, {phiE, phiE}, [3, 2]));
        f11E = 2 * f20E;
        f21E = ttv(dT2dT2, {w20E + w11E, phiE}, [3, 2]) + ...
               3 * ttv(dT3, {phiE, phiE, phiE}, [4, 3, 2]);
        f21E = double(f21E);
        
        % Compute sensitivity
        dgamma(ii, :) = 1/(2*omega) * phiE.' * f21E + ...
                        lambda0E.' * (dK - omega^2 * dM) * phiE + ...
                        lambda1 * phiE.' * dM * phiE + ...
                        lambda20E.' * (dL20 * w20E + f20E) + ...
                        lambda11E.' * (dL11 * w11E + f11E);
    end
end
