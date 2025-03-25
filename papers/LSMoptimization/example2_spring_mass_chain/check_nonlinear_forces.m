% Compute the LSM of a spring-mass chain
clear; clc; close all;
addpath(genpath(fullfile('..', 'LSM-opt-functions')))

%% Settings
% Define the type of the system
sysArgs.systemType = @SystemSpringMassChain;

% Define the boundary conditions: 0 to fix the mas, 1 to let it free.
sysArgs.bc = [0, ones(1, 3), 0];
nMasses = sum(sysArgs.bc);
nSprings = nMasses + 1;

% Define the system parameters
% m = 3;
% k = 30;
% k2 = 3;
% k3 = 3;
% x0u = [m, k, k2, k3];
% x = [x0u(1) * ones(1, nMasses), ...
%      x0u(2) * ones(1, nSprings), ...
%      x0u(3) * ones(1, nSprings), ...
%      x0u(4) * ones(1, nSprings)];
m  = rand(1, nMasses);
k  = rand(1, nSprings);
k2 = rand(1, nSprings);
k3 = rand(1, nSprings);
x = [m, k, k2, k3];

% Build system
sys = feval(sysArgs.systemType, x, sysArgs);

%% Validate nonlinear forces

% Random displacement
N = length(sysArgs.bc);
n = N - 2;
w = rand(N, 1);
w(1) = 0; w(end) = 0;
z = [w(2:end-1, 1); zeros(n, 1)];

% Compute the quadratic nonlinearity
f2 = zeros(size(w));
f2(1) = -k2(1) * (w(2) - w(1))^2;
for i = 2:length(w)-1
    f2(i) = k2(i-1) * (w(i) - w(i-1))^2 - k2(i) * (w(i+1) - w(i))^2;
end
f2(end) = k2(end) * (w(end) - w(end-1))^2;

% Compute the cubic nonlinearity
f3 = zeros(size(w));
f3(1) = -k3(1) * (w(2) - w(1))^3;
for i = 2:length(w)-1
    f3(i) = k3(i-1) * (w(i) - w(i-1))^3 - k3(i) * (w(i+1) - w(i))^3;
end
f3(end) = k3(end) * (w(end) - w(end-1))^3;

% Compute the tensor contraction
f2t = sys.f{2} * kronexp(z, 2);
f3t = sys.f{3} * kronexp(z, 3);

% Errors
norm(f2(2:end-1) - f2t)
norm(f3(2:end-1) - f3t)
