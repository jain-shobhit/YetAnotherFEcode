% HB_sensitivity
%
% synthax: SENS = HB_sensitivity(X, systemNom, derFunHandle, m, method, ...
%                                 approxOrder, Ns, varargin)
% DESCRIPTION:
% HB_sensitivity function returns with the sensitivity analysis of HB
% coefficients with respect to some set of parameters. Sensitivity
% technique implemented here is described in [2]
% 
% Moreover it must hold:
% (1)E.o.m. in time domain for the studied problem must be in the form:
%    Mx_dd + Dx_d + K(p)x + f_nl(x,p) + f_ex = 0. 
% (2)Internal forces must be nonlinear only in x, (not in x_d) 
% (3)parameter dependent terms are only internal stiffness forces 
%    (M and D independent from parameter vector p)
% (These constraints come from function 'HBderivatives' used in this code)
% 
% INPUTS:
% (1) X:            column vector containing linearization point. It is
%                   often the point of solution of the nominal problem.
%                   X contains HB coefficients, and last element contains
%                   angular frequency omega (order of elements as in NlVib
%                   [1])
% (2) systemNom:    structure array containing fields 'K','M','D','n',
%                   identifying the stiffness, mass and damping matrices of
%                   the nominal system, and number of dof.
% (3) derFunHandle: derFunHandle =  @(q)fderivatives(q,other_arguments).
%                   fderivatives(q) evaluates internal forces derivatives 
%                   as a function of the vector of generalized displacements
%                   q. Output of fderivatives(q) is a struct that must 
%                   contain the following fieds: dfdp, dfdq, dfdq2, dfdp2,
%                   dfdqdp (i.e. derivatives of f wrt q and parameters p).
% (4) m:            size of parameter vector p
% (5) method:       possible options are 'normal' (1) or 'omegaconst' (2).
%                   case (1) allows to perform sensitivity analysis with 
%                   the normal method (see ref *). Unless specified in 
%                   varargin{1} weight on \delta\omega = 0, while the other
%                   coefficients have the same weight. 
%                   case (2) is used for sensitivity at constant omega.
% (6) approxOrder:  possible options are '1' (1) or '2' (2). case (1) only
%                   first order sensitivities are computed, case (2) first
%                   and second order sensitivities are computed.
% (7) Ns:           number of samples used in the AFT algorithm to compute
%                   derivatives of the residual.   
% (8) varargin{1}:  if scalar varargin{1} specifies the value of the 
%                   relative weight of omega with respect to HB 
%                   coefficients, that weight the same (1). If you want to 
%                   give customize weight to HB coefficients and omega fill 
%                   varargin{1} with diagonal matrix of scaling factors 
%                   (same size of X).
%
% OUTPUTS:
% (1) SENS:         structure array with fields 'S1' (if approxOrder == 1) 
%                   or 'S1' and 'S2' (if approxOrder == 2). SENS.S1 
%                   contains first order sensitivities, SENS.S2 contains      
%                   second order sensitivities. 
%
% FURTHER IMPLEMENTATIONS:
% all systems share same matrix to be inverted. Find a way to substitute 
% inversion made with '\' with the inverse computed only once without
% loosing accuracy.
%
% REFERENCES:
% [1]: NLvib Version 1.3, Copyright (C) 2020 Malte Krack, Johann Gross
% [2]: Saccani... 
%
% Author: Alexander Saccani, Msc in mechanical engineering
% University: Politecnico di Milano
% Created:  09/2021

function SENS = HB_sensitivity(X, systemNom, derFunHandle, m, method, ...
                                approxOrder, Ns, varargin)

% set input
if nargin < 8
    W = eye(length(X));
    W(end,end) = 0;
elseif isnumeric(varargin{1}) %check input
    if size(varargin{1}) == 1
        scalingOmega = varargin{1};
        W = eye(length(X));
        W(end,end) = scalingOmega;
    else
        W = varargin{1};
    end
else
    error('input varargin{1} must be double')
end



% *************************************************************************
% Evaluate partial derivatives of the residual
% *************************************************************************

n = systemNom.n;                % number of dofs
H = ((length(X)-1)/n -1 )/2;    % number of harmonics

% cosine and sine indeces
IC = n+repmat(1:n,1,H)+n*kron(0:2:2*(H-1),ones(1,n));
IS = IC+n; 

% reorganize HB coefficients in 3D matrix
% (Xt)_{ijk}:
%   i -> dof number
%   j -> harmonic number
%   k -> 1=A, 2=B
Xt = zeros(n, H+1, 2); 
Xt(:, 1, 1) = X(1:n); % static term (0th harmonic)
for h = 1:H 
    Xt(:, h+1, 1) = X( IC((h-1)*n+1 : (h-1)*n+n) );
    Xt(:, h+1, 2) = X( IS((h-1)*n+1 : (h-1)*n+n) );
end

Om = X(end); %omega

derivatives = HB_derivatives(Xt, Om, systemNom, derFunHandle, m, Ns);

% rename derivatives
dRdP   = derivatives.dRdP;
dRdX   = derivatives.dRdX;
dRdXdX = derivatives.dRdX2;
dRdPdP = derivatives.dRdP2;
dRdXdP = derivatives.dRdXdP;
dRdPdX = permute(dRdXdP,[1,3,2]);

dimX = length(X);       % dimension of HB coeff vector
dimP = size(dRdP,2);	% dimension of parameter vector
if dimP ~= m
    error(['decleared number of paremeters differs from number of ',...
        'parameters in user provided derivatives']);
end



% *************************************************************************
% Perform sensitivity analysis
% *************************************************************************

switch(lower(method))
    case('normal')
        
        % FIND TANGENT DIRECTION __________________________________________
        C = dRdX;
        % find tangent vector
        C_Q = C(1:end,1:end-1);
        C_om = C(:,end);
        t = [-C_Q\C_om;1]; % tangent vector
        % apply scaling to tangent vector
        t_til = t.*diag(W.^2); 

        % COMPUTE 1ST ORDER SENSITIVITY ___________________________________
        C_til = [dRdX; t_til'];
        D_til = [dRdP; zeros(1,size(dRdP,2))]; %known term
        S1 = -C_til\D_til; %first order sensitivity
        
        % SECOND ORDER SENSITIVITIES ______________________________________
        if approxOrder == 2
            % CONSTRUCT TENSORS FOR SECOND ORDER SENSITIVITY
            % convert to tensor
            S1t = tensor(S1);
            
            if m == 1
                G = squeeze(ttt(tensor(dRdXdX),S1t,3,1)) + tensor(dRdXdP);
                F = squeeze(ttt(tensor(dRdPdX),S1t,3,1)) + squeeze(tensor(dRdPdP));
                E = ttt(G,S1t,2,1);
            else
                G = ttt(tensor(dRdXdX),S1t,3,1) + tensor(dRdXdP);
                F = ttt(tensor(dRdPdX),S1t,3,1) + tensor(dRdPdP);
                E = ttt(G,S1t,2,1);
            end
            
            % convert back to double
            F = double(F);
            E = double(E);
            
            % compute derivative of tangent vector
            L = double(ttv(G,t,2));
            dRdQ = dRdX(:,1:end-1);
            dtdP = eye(size(W,1)).*diag(W.^2)*[-dRdQ\L;zeros(1,m)]; %dtdOm = 0
            M = -ttt(tensor(dtdP),S1t,1,1);
            M = double(M);

            % COMPUTE SECOND ORDER SENSITIVITIES___________________________

            %initialize second order sensitivity
            S2 = zeros(dimX, dimP, dimP);

            %solve second order sensitivity systems of equations
            kEnd = dimP;
            kStart = 1;
            
            % for schwartz theorem, no need to compute all 2nd ord. deriv.
            for jj = 1:dimP
                for kk = kStart:kEnd

                    C_til2 = [C; t_til'];
                    D_til2 = [F(:,jj,kk) + E(:,kk,jj); M(jj,kk)]; 

                    S2jjkk = -C_til2\D_til2;

                    S2(:,jj,kk) = S2jjkk;
                    S2(:,kk,jj) = S2jjkk;

                end
                kStart = kStart + 1;
            end
        end
      
    case('omegaconst')
    
        dRdQ = dRdX(:,1:end-1);
        dRdQdQ = dRdXdX(:,1:end-1,1:end-1);
        dRdPdQ = dRdPdX(:,:,1:end-1);
        dRdQdP = permute(dRdPdQ,[1 3 2]);
        
        % COMPUTE 1ST ORDER SENSITIVITY____________________________________
        C = dRdQ;
        D = dRdP; %known term
        S1 = -C\D;
        S1 = [S1;zeros(m,1)]; % add sensitiviy of omega
        
        % SECOND ORDER SENSITIVITIES_______________________________________
        if approxOrder == 2
            
            S1t = tensor(S1(1:end-1,:));
            
            if m == 1 
                G = squeeze(ttt(tensor(dRdQdQ),S1t,3,1)) + tensor(dRdQdP);
                F = squeeze(ttt(tensor(dRdPdQ),S1t,3,1)) + squeeze(tensor(dRdPdP));
                E = ttt(G,S1t,2,1);
            else
                G = ttt(tensor(dRdQdQ),S1t,3,1) + tensor(dRdQdP);
                F = ttt(tensor(dRdPdQ),S1t,3,1) + tensor(dRdPdP);
                E = ttt(G,S1t,2,1);
            end
            
            %convert back to double
            F = double(F);
            E = double(E);

            % COMPUTE SECOND ORDER SENSITIVITIES___________________________
            dimX = length(X);
            dimQ = length(X)-1; %dimension of HB coeff vector
            dimP = size(dRdP,2); %dimension of parameter vector

            %initialize second order sensitivity
            S2 = zeros(dimQ,dimP,dimP);

            %solve second order sensitivity systems of equations
            kEnd = dimP;
            kStart = 1;

            for jj = 1:dimP
                for kk = kStart:kEnd

                    D2 = F(:,jj,kk) + E(:,kk,jj); 

                    S2jjkk = -C\D2;

                    S2(:,jj,kk) = S2jjkk;
                    S2(:,kk,jj) = S2jjkk;

                end
                kStart = kStart + 1;
            end

            S2tmp = S2;
            S2 = zeros(dimX,dimP,dimP);
            S2(1:end-1,:,:) = S2tmp;            
            
        end
    
    otherwise
        error(['Non valid method for sensitivity computation. Possible ',...
            'options are: normal (1), omegaconst (2)']);
end



% *************************************************************************
% Generate output
% *************************************************************************

% output generation
SENS.Sens.S1 = S1;

if approxOrder == 2
    SENS.Sens.S2 = S2;
end

end