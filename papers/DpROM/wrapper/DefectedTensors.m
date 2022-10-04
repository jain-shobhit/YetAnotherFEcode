% DefectedTensors
%
% Synthax:
% [Q2, Q3, Q4, Q3t, Q4t, M] = DefectedTensors(tensors, xi)
%
% Description: this function returns the second, third and fourth order
% reduced stiffness tensors evaluated at xi (see Eqs. (28) and (33) in the
% paper).
% INPUTS
%   - tensors: the output of the "reduced_tensors_DpROM.m" function
%   - xi: column vector containing the amplitudes of the defects in the
%     selected basis U (we'll have then: d=size(U,2)=length(xi))
% OUTPUTS
%   - Q2, Q3 and Q4: the reduced order stiffness tensors
%   - Q3t, Q4t: tensors required to compute the reduced tangent stiffness
%     matrix (see tensors_KF.m function)
%   - M: parametrized mass matrix, computed using an approximated
%     integration over the defected volume. The defected mass matrix can be
%     obtained as:
%           Md = M{1} + M{2}*xi(1) + M{3}*xi(2) + ... M{d+1}*xi(d)
%
% Reference: J. Marconi, P. Tiso, D.E. Quadrelli & F. Braghin, "A higher 
% order parametric nonlinear reduced order model for imperfect structures 
% using Neumann expansion", Nonlinear Dynamics, 2021.
%
% Created: 14 May 2021
% Last modified: 21 October 2021
% Author: Jacopo Marconi, Politecnico di Milano

function [Q2, Q3, Q4, Q3t, Q4t, M] = DefectedTensors(tensors, xi)

nd = length(xi);

% tensors computed on the nominal volume
Q2n = tensors.Q2n{1};
Q3n = tensors.Q3n{1};
Q4n = tensors.Q4n{1};
Q3d = tensors.Q3d{1};
Q4d = tensors.Q4d{1};
Q5d = tensors.Q5d{1};
Q4dd = tensors.Q4dd{1};
Q5dd = tensors.Q5dd{1};
Q6dd = tensors.Q6dd{1};
M = tensors.M{1};

% apply volume correction to integrate over the defected volume
if tensors.volume == 1
    for dd = 2 : nd+1
        Q2n  = Q2n  + tensors.Q2n{dd}  * xi(dd-1);
        Q3n  = Q3n  + tensors.Q3n{dd}  * xi(dd-1);
        Q4n  = Q4n  + tensors.Q4n{dd}  * xi(dd-1);
        Q3d  = Q3d  + tensors.Q3d{dd}  * xi(dd-1);
        Q4d  = Q4d  + tensors.Q4d{dd}  * xi(dd-1);
        Q5d  = Q5d  + tensors.Q5d{dd}  * xi(dd-1);
        Q4dd = Q4dd + tensors.Q4dd{dd} * xi(dd-1);
        Q5dd = Q5dd + tensors.Q5dd{dd} * xi(dd-1);
        Q6dd = Q6dd + tensors.Q6dd{dd} * xi(dd-1);
        M    = M    + tensors.M{dd}    * xi(dd-1);
    end
end

if isscalar(xi)
    % slightly different syntax if xi is a scalar
    Q3d  = Q3d*xi;
    Q4d  = Q4d*xi;
    Q5d  = Q5d*xi;
    Q4dd = Q4dd*xi^2;
    Q5dd = Q5dd*xi^2;
    Q6dd = Q6dd*xi^2;
else
    Q3d  = ttv(Q3d,xi,3);
    Q4d  = ttv(Q4d,xi,4);
    Q5d  = ttv(Q5d,xi,5);
    Q4dd = ttv(ttv(Q4dd,xi,4),xi,3);
    Q5dd = ttv(ttv(Q5dd,xi,5),xi,4);
    Q6dd = ttv(ttv(Q6dd,xi,6),xi,5);
end

Q2 = Q2n + double(Q3d + Q4dd);
Q3 = Q3n + double(Q4d + Q5dd);
Q4 = Q4n + double(Q5d + Q6dd);

% for the tangent stiffness matrix
Q3t = Q3 + permute(Q3, [1 3 2]); 
Q4t = Q4 + permute(Q4, [1 3 2 4]) + permute(Q4, [1 4 2 3]);
