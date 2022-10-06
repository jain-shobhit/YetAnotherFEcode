% DESCRIPTION
% this function projects to full space the reduced coordinates
%
% INPUTS:
% (1) V:             reduction base of pROM
% (2) frc_red_coord: Struct array with fields, Q0, Qre, Qim,
%                    omega in which the reduced HB coefficients are saved.
%                    Structure of frc_red_coord is the same of output of
%                    'nlvib_decode.m'.
%
% OUTPUTS:
% (1) full_disp:     HB coefficients of full coordinates. 
%
% Author: Alexander Saccani, Msc in mechanical engineering
% University: Politecnico di Milano
% Created:  09/2021

function full_disp = get_full_disp(V, frc_red_coord)

N_samples = length(frc_red_coord.omega);
N_dof_full = size(V, 1);
H = size(frc_red_coord.Qre, 3);

%initialize output
full_disp.omega = frc_red_coord.omega;
full_disp.Qre = zeros(N_dof_full, N_samples, H);
full_disp.Qim = zeros(N_dof_full, N_samples, H);
full_disp.Q0  = zeros(N_dof_full, N_samples);

%fill matrices
full_disp.Q0 = V*frc_red_coord.Q0;

for n_H = 1:H
    full_disp.Qre(:,:,n_H) = V * frc_red_coord.Qre(:, :, n_H);
    full_disp.Qim(:,:,n_H) = V * frc_red_coord.Qim(:, :, n_H);
end
    
end

