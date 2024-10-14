function [c, ceq, dc, dceq] = template_nonlcon(x, sysArgs, lsmArgs, optArgs, doPrint)
% Nonlinear constraints function for optimization. This function computes
% the equality and inequality constraints and their sensitivities.
% Inputs:
%   x: the design variables.
%   sysArgs: the arguments structure.
%   lsmArgs: the LSM arguments structure.
%   optArgs: the optimization arguments structure.
%   doPrint: a flag to print the computation times. Default is false.
% Outputs:
%   c: the inequality constraints.
%   ceq: the equality constraints.
%   dc: the inequality constraints sensitivities.
%   dceq: the equality constraints sensitivities.

if nargin == 4
    doPrint = false;
end

% Build system
t0=tic;
sys = feval(sysArgs.systemType, x, sysArgs);
if doPrint
    fprintf(' System built in %.1f s\n', toc(t0))
end

% Modal analysis
sys = sys.modal_analysis();

% Identify master mode
if isfield(sysArgs, 'modeIndex')
    % Assign master mode
    modeIndex = sysArgs.modeIndex;
elseif isfield(sysArgs, 'modeReference')
    % Apply modal assurance criterion
    modeIndex = sys.modal_assurance_criterion(sys.phi, sysArgs.modeReference);
else
    % Default master mode
    modeIndex = 1;
end

% Assign master mode
sys.omega0 = sys.omega(modeIndex);
sys.phi0 = sys.phi(:, modeIndex);
sys = sys.state_space_modal_analysis();

% Compute LSM
t0=tic;
lsm = LSM(sys, lsmArgs.order);
lsm = lsm.compute_manifold();
if doPrint
    fprintf(' LSM computed in %.1f s\n', toc(t0))
end

% Compute LSM sensitivities
t0=tic;
sens = Sensitivity(sys, lsm);
sens = sens.sensitivity_lsm();
if doPrint
    fprintf(' Sensitivities computed in %.1f s\n', toc(t0))
end

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
