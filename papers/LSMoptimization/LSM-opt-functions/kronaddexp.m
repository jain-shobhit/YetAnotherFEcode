function out = kronaddexp(in, n)
% Performs the Kronecker tensor sum on variable "in" "n" times
% Inputs:
%   in: the input matrix.
%   n: the number of times the Kronecker tensor sum is performed.
% Outputs:
%   out: the output matrix.

out = in;
for ii = 2:n
    out = kronadd(out, in);
end
