function out = kronn(ins)
% Performs the Kronecker tensor product between all the elements in the cell array "ins"
% Inputs:
%   ins: cell array, input matrices or vectors.
% Outputs:
%   out: the output matrix or vector.

out = ins{1};
for ii = 2:length(ins)
    out = kron(out, ins{ii});
end
