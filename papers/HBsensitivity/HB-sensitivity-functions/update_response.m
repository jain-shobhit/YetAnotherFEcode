% DESCRIPTION:                                                      
% update_response allows to obtain the varied HB coefficients for a 
% perturbation in the parameter vector dP from its nominal value. Update is
% sensitivity based (linear approximation or quadratic approximation
% according if only first or first and second order sensitivities are
% provided). 
%
% INPUTS:
% (1) Xnom:    double array containinig samples of nominal HB coefficients. 
%              It is the output of NlVib (see ref*). Nominal solution is
%              computed for dP = 0.
% (2) Sens:    Structure array with fields 'S1' (1) or 'S1' and 'S2'(2). In
%              case (1) update_response will provide you with the perturbed
%              solution obtained with linear sensitivities while in case
%              (2) the function will provide both linear and quadratic
%              approximations (1st and 2nd order Taylor expansions).
% (3) dP:      array containing parameter variation vector.
%
% OUTPUTS:
% (1) Xup:     Structure array containing sensitivity based perturbed HB 
%              coefficents for parameter variation dP. Fields are 'lin' (1)
%              or 'lin' and 'quad' (2) depending on input variable 'Sens'.
%
% Author: Alexander Saccani, Msc in mechanical engineering
% University: Politecnico di Milano
% Created:  09/2021

function Xup = update_response(Xnom,Sens,dP)


% initialize output
n_samples = size(Xnom,2);       % number of samples
dimX = size(Xnom,1);            % size of X
m = length(dP);                 % length of parameter vector
Xup1 = zeros(dimX,n_samples);
var1 = zeros(dimX,n_samples);

% first order update
for ii = 1:n_samples
    S1 = Sens(ii).S1;
    var1(:,ii) = S1*dP;
    Xup1(:,ii) = Xnom(:,ii) + var1(:,ii);
end

% second order update
flag2 = 0;
if isfield(Sens,'S2')
    flag2 = 1;
    Xup2 = zeros(dimX,n_samples); %inizialize output
    for ii = 1:n_samples        
        S2 = Sens(ii).S2;
        if m == 1
            var2 = 0.5*S2*dP^2;
        else
            var2 = 0.5*ttv(tensor(S2),{dP,dP},[2,3]);
            var2 = double(var2);
        end        
        Xup2(:,ii) = Xnom(:,ii) + var1(:,ii) + var2;         
    end    
end

% generate output
Xup.lin = Xup1;
if flag2
    Xup.quad = Xup2;
end

end

