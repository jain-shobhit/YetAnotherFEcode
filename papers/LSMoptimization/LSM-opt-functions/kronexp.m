function out = kronexp(in, n)
% Performs the Kronecker tensor product on variable "in" "n" times
% Inputs:
%   in: the input matrix or vector.
%   n: the number of times the Kronecker product is performed.
% Outputs:
%   out: the output matrix or vector.

out = in;
for ii = 2:n
    out = kron(out, in);
end
