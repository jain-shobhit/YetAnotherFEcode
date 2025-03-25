% Check the nonlinear forces of a uniform spring-mass chain
clear; clc; close all;
addpath(genpath(fullfile('..', 'LSM-opt-functions')))

%% Settings
% Define the type of the system
sysArgs.systemType = @SystemSpringMassChainUniform;

% Define the boundary conditions: 0 to fix the mas, 1 to let it free.
sysArgs.bc = [0, ones(1, 5), 0];

% Define the system parameters
m = 3;
k = 30;
k2 = 3;
k3 = 3;
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
f2(1) = -k2 * (w(2) - w(1))^2;
for i = 2:length(w)-1
    f2(i) = k2 * (w(i) - w(i-1))^2 - k2 * (w(i+1) - w(i))^2;
end
f2(end) = k2 * (w(end) - w(end-1))^2;

% Compute the cubic nonlinearity
f3 = zeros(size(w));
f3(1) = -k3 * (w(2) - w(1))^3;
for i = 2:length(w)-1
    f3(i) = k3 * (w(i) - w(i-1))^3 - k3 * (w(i+1) - w(i))^3;
end
f3(end) = k3 * (w(end) - w(end-1))^3;

% Compute the tensor contraction
f2t = sys.f{2} * kronexp(z, 2);
f3t = sys.f{3} * kronexp(z, 3);

% Errors
norm(f2(2:end-1) - f2t)
norm(f3(2:end-1) - f3t)
