function rms = RMS(frc_curve,harm)
% DESCRIPTION
% RMS function to get root mean squared value of frc for given harmonics.
%  
% INPUTS:
% (1) frc_curve: struct array containing frc curve. Array must contain the
%                following fields: 'Q0', 'Qre', 'Qim', 'omega'. frc_curve
%                could be output of 'nlvib_decode.m'
% (2) harm:      vector array containing all harmonics that you want to
%                include in the rms
% OUTPUTS: 
% (1) rms: double array containing root mean squared value of frc modulus 
%          for each dof (for specified harmonics). 
%          rms(ii,kk) -> rms(dof(ii),frequency sample kk)
%
% Author: Alexander Saccani, Msc in mechanical engineering
% University: Politecnico di Milano
% Created:  09/2021

Nsamples = length(frc_curve.omega);
Ndofs = size(frc_curve.Q0,1);

% initialize output
rms = zeros(Ndofs,Nsamples);

% contribution of static term
flag = 0;
if ~isempty(find(harm == 0))
    
    flag = 1;
    harm(find(harm == 0)) = [];
    
    for n_dof = 1:Ndofs  
        
        Q0 = frc_curve.Q0;
        Q0 = Q0(n_dof,:);
        rms(n_dof,:) = Q0.^2;
        
    end
    
    if isempty(harm)         
        rms = sqrt(rms);
        return
    end
    
end

% add contribution of other harmonics
for n_dof = 1:Ndofs
    
   for ii = 1:length(harm)
       
       n_harm = harm(ii);
       A = frc_curve.Qre(:,:,n_harm);
       B = frc_curve.Qim(:,:,n_harm);
       Q = A + 1i*B;
       Q = Q(n_dof,:);
       mod_harm = abs(Q);
       rms(n_dof,:) = rms(n_dof,:) + mod_harm.^2;
       
   end
   
   if flag == 1
    rms(n_dof,:) = sqrt(rms(n_dof,:)/( length(harm)+1 ));
   else
    rms(n_dof,:) = sqrt(rms(n_dof,:)/length(harm));
   end
   
end

end

