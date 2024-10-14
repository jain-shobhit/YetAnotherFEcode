function err = error_measure(lsm, lsmArgs, optArgs)
% Calculate the error measure [Eq. (27)] from the invariance equation.
% Inputs:
%   lsm: LSM object.
%   lsmArgs: LSM arguments.
%   optArgs: optimization arguments.
% Outputs:
%   err: error measure.

% Maximum of the target displacement and corresponding rho
zTarget = optArgs.zTarget;
rho = lsm.solve_rho(zTarget(end), lsmArgs.iDof);

% Number of theta values
nTheta = 30;
theta = linspace(0, 2*pi, 30);

% Loop over theta values
err = 0;
for ii = 1:nTheta
    [resNorm, norms] = invariance_equation_residual(lsm, rho, theta(ii), lsm.order);
    err = err + resNorm / nTheta / max(norms);
end
