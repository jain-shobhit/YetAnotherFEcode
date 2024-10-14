function [resNorm, norms] = invariance_equation_residual(lsm, rho, theta, order)
% Computes the residual of the invariance equation
% res = B * DpW(p) * R(p) - A * W(p) - F(W(p))
% Inputs:
%   lsm: the Lyapunov Subcenter Manifold object.
%   rho: the radial coordinate.
%   theta: the angular coordinate.
%   order: the order of the expansion.
% Outputs:
%   resNorm: the residual's norm.
%   norms: the norms of the terms of the residual.

% Extract variables
W = lsm.W;
R = lsm.R;
F = lsm.sys.F;
A = lsm.sys.A;
B = lsm.sys.B;

% Check inputs
if nargin < 4
    order = length(W);
end

% Point in the reduced space
p = rho * exp([1i*theta; -1i*theta]);

% Evaluate W(p) and R(p)
pExp = p;
Rp = R{1} * pExp;
Wp = W{1} * pExp;
for ii = 2:order
    pExp = kron(pExp, p);
    Rp = Rp + R{ii} * pExp;
    Wp = Wp + W{ii} * pExp;
end

% F(W(p))
WpExp = Wp;
FWp = F{1} * WpExp;
for ii = 2 : length(F)
    WpExp = kron(WpExp, Wp);
    FWp = FWp + F{ii} * WpExp;
end

% B*DpW(p)*R(p)
DpW = 0;
for ii = 1 : order
    DpW = DpW + W{ii} * differential_kronexp(p, ii);
end
BDWR = B*DpW*Rp;

% Residual's norm
resNorm = norm(BDWR - A*Wp - FWp);
norms = [norm(BDWR), norm(A*Wp), norm(FWp)];
