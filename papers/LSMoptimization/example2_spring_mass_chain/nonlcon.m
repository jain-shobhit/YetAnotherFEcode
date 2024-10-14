function [c, ceq, dc, dceq] = nonlcon(x, sysArgs, lsmArgs, optArgs)
% Nonlinear constraints function for optimization. This function computes
% the equality and inequality constraints and their sensitivities.
% Inputs:
%   x: the design variables.
%   sysArgs: the arguments structure.
%   lsmArgs: the LSM arguments structure.
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

% Default master mode
modeIndex = sysArgs.modeIndex;

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
