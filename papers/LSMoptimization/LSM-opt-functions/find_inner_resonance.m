function [lIndices, jIndices, permIndices, lambdaFull] = find_inner_resonance(order)
% Find the indicesof the inner resonances at a given order.
% Inputs:
%   order: scalar, the order of expansion.
% Outputs:
%   lIndices: vector, the indices of the inner resonances.
%   jIndices: vector, the indices of the inner resonances.
%   Delta: matrix, the indices of the inner resonances.
%   lambdaFull: vector, the eigenvalues of the full system at the given order

% Adimensional eigenvalues
lambda = [1, -1];

% Number of modes
M = 2;

% Build all the permutations for the given order
v = (1:M).';
permIndices = disp_with_rep(v, order);

% Compute the eigenvalues of the full system at the given order
lambdaFull = zeros(M^order, 1);
for ii = 1:size(permIndices, 1)
    lambdaFull(ii) = sum(lambda(permIndices(ii, :)));
end

% Identify inner resonances
lambda1 = lambdaFull == 1;
lambda2 = lambdaFull == -1;

% Inner resonance indices
temp = (1:(M^order)).';
lIndices = [temp(lambda1); temp(lambda2)];
jIndices = [1*ones(sum(lambda1), 1); 2*ones(sum(lambda2), 1)];
