% tensors_KF
%
% Synthax:
% [Kt,fi] = tensors_KF(Q2, Q3, Q4, Q3t, Q4t, q)
%
% Description: computes reduced tangent stiffness matrix and reduced
% internal forces using the tensorial formulation.
% INPUTS
%   - Q2, Q3, Q4: reduced stiffness tensors of order 2, 3 and 4
%   - Q3t, Q4t: tensors required to compute the tangent stiffness matrix.
%     They can be obtained from Q3 and Q4 as:
%       Q3t = Q3 + permute(Q3, [1 3 2]); 
%       Q4t = Q4 + permute(Q4, [1 3 2 4]) + permute(Q4, [1 4 2 3]);
%     Since, however, Q3t and Q4t are constant, it is more efficient to
%     compute them outside this function once and for all.
%   - q: reduced coordinates vector
% OUTPUTS
%   - Kt: reduced tangent stiffness matrix
%   - fi: reduced internal forces
%
% Created: 14 April 2021
% Author: Jacopo Marconi, Politecnico di Milano

function [Kt,fi] = tensors_KF(Q2,Q3,Q4,Q3t,Q4t,q)

fi = Q2*q + ttsv(Q3,q,-1) + ttsv(Q4,q,-1);
Kt = Q2 + ttsv(Q3t, q, -2) + ttsv(Q4t,q,-2);