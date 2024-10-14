function d = disp_with_rep(v, m)
% Generate all the m dispositions of the elements in the vector v.
% Inputs:
%   v: vector, the elements to be disposed (n x 1)
%   m: scalar, the number of dispositions.
% Outputs:
%   d: matrix, the dispositions.

% Number of elements
n = length(v);

% Number of dispositions
nDisp = n^m;

% Compute dispositions
d = zeros(nDisp, m);
for ii = 1:m
    d(:, ii) = kron(kron(ones(n^(ii - 1), 1), v), ones(n^(m - ii), 1));
end
