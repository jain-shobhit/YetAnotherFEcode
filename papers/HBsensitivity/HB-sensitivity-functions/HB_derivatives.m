% DESCRIPTION:                                                      
% HB_derivatives allows to numerically evaluate the first and second order
% derivatives of the HB residual with respect to HB coefficients and to a
% specified parameter vector p at specified solution point Xt. Derivatives
% of the residual are computed from time domain snapshots with the
% alternate Fourier transform (A.F.T.) algorithm.
% 
% Moreover it must hold:
% (1)E.o.m. in time domain for the studied problem must be in the form:
%    Mq_dd + Dq_d + K(p)q + f_nl(q,p) + f_ex = 0. 
% (2)Internal forces must be nonlinear only in q, (not in q_d) 
% (3) parameter dependent terms are nonlinear internal stiffness forces. 
%  M, D and K can be only linearly dependent on p. Anyway it is still 
%  possible to assume any particular dependence of K on p (not just
%  linear) by considering it as part of nonlinear force vector.
%  (put dKdp = 0, and consider in dfdp, dfdp2... the derivatives of linear
%  forces in displacements.
% 
% INPUTS: 
% (1)Xt:           multidimensional matrix containing solution point 
%                  in HB coefficients space. Xt{ijk} i-> dof, j-> harm-1,
%                  if k == 1 -> A, if k == 2 -> B. 
% (2)Om:           scalar array containing angular frequency value.
% (3)system:       structure array containing fields 'M','K','D' that
%                  are the mass, stiffness and damping matrix of the system
% (4)JacFunHandle: JacFunHandle =  @(q)fderivatives(q,other_arguments).
%                  fderivatives(q) evaluates internal forces derivatives 
%                  as a function of the vector of generalized displacements
%                  q. Output of fderivatives(q) is a struct that must 
%                  contain the following fieds: dfdp, dfdq, dfdq2,dfdp2,
%                  dfdqdp. If Mass, damping and stiffness matrices depend
%                  on parameters then function must contain fields dMdp,
%                  dKdp, dDdp. If dMdp, dKdp, dDdp are not provided fileds,
%                  then they are consider null. It is compulsary that K,M,D
%                  are just linearly dependent on p. With a trick it is
%                  also possible to consider K dependt on p with any
%                  function law (not necessarely linear) (see point 3 of
%                  function description).
%     
% (5)m:            lenght of parameter vector p
% (6)Ns:           number of samples used in the AFT algorithm to compute
%                  derivatives of the residual. 
%
% OUTPUTS:
% (1) HBderivatives: structure array with the following fields:
%                     'dRdX', 'dRdP', 'dRdX2', 'dRdP2', 'dRdXdP'
%                   order of elements in dimension X corresponds to order
%                   of HB coefficients and omega implemented in NlVib tool
%                   [1].
%
% POSSIBLE FUTURE IMPLEMENTATIONS:
% (1) according to Schwarz's theorem, second order derivatives are
% symmetric. No need to compute them two times. Reduce computational
% effort!
% (2) implement this function in such a way that nonlinear forces can
% depend also on velocities f_nl = f_nl(x,x_d,p)
%
%
% REFERENCES:
% [1]: NLvib Version 1.3, Copyright (C) 2020 Malte Krack, Johann Gross
% 
% Author: Alexander Saccani, Msc in mechanical engineering
% University: Politecnico di Milano
% Created:  09/2021

function HBderivatives = HB_derivatives(Xt, Om, system, JacFunHandle, m, Ns)
%% evaluation of derivatives in time domain (X->q(t)->dfdq(t),...   
n = size(Xt,1); %number of dof
H = size(Xt,2)-1; %remove static term

M = system.M;
K = system.K;
D = system.D;

% evaluate displacements in time domain
T = 2*pi/Om;
dt = T/Ns; %time sample for fft
t = 0:dt:T-dt; %time snapshots for fft

Qce = zeros(Ns,n);
for h = 1:H+1
    Qce(h,:) = 0.5*Xt(:,h,1).' -0.5*1i*Xt(:,h,2).';
end
q = 2*real(Ns*ifft(Qce).'); %columns indicate time instant, rows dof. q contains displacements evaluations during the period

%extract derivatives values of forces in time domain
derivatives = extend_der_in_time(q,JacFunHandle);
dMdKdD = feval(JacFunHandle,q(:,1)); %to find derivatives of K,M,D that do not depend on q.

dfdq = derivatives.dfdq; %only nl forces
dfdp = derivatives.dfdp; %linear + nl forces
dfdq2 = derivatives.dfdq2; %only nl forces (linear forces der is null)
dfdqdp = derivatives.dfdqdp; %linear + nl forces
dfdp2 = derivatives.dfdp2; %linear + nl forces

%extract derivatives of mass, stiffness and damping matrices
if isfield(dMdKdD,'dMdp')
    dMdp = tensor(dMdKdD.dMdp);
else
    dMdp = tensor(zeros(n,n,m)); %improve efficiency!!
end

if isfield(dMdKdD,'dKdp')
    dKdp = tensor(dMdKdD.dKdp);
else
    dKdp = tensor(zeros(n,n,m));
end

if isfield(dMdKdD,'dDdp')
    dDdp = tensor(dMdKdD.dDdp);
else
    dDdp = tensor(zeros(n,n,m));
end


%% dRdXt                                                            
dRdXt = zeros(n,H+1,2,n,H+1,2); %need to consider also static terms (harm_number = H+1)
dRdXt_nl = dRdXt;

%fill with linear contributions
for J = 1:H+1
    
    %derive cosine terms of ft of residual
    dRdXt(:,J,1,:,J,1) = -M*Om^2*(J-1)^2 + K; %derive w.r.t. A
    dRdXt(:,J,1,:,J,2) = D*Om*(J-1); %derive w.r.t. B
    
    %derive sine terms of ft of residual
    dRdXt(:,J,2,:,J,1) = -D*Om*(J-1);
    dRdXt(:,J,2,:,J,2) = -M*Om^2*(J-1)^2 + K;

end

%compute jacobian of nonlinear terms with AFT
Mcosine = cos( (0:1:H)'*Om*t );
Msine = sin( (0:1:H)'*Om*t );

for I = 1:n
    
    for L = 1:n
        dfI_dqL = squeeze(dfdq(I,L,:))'; %vector of function samples in time domain, last dimension is time sample
        
        for N = 1:H+1 %for every harmonic except the first one
            
            fft_1 = fft(dfI_dqL.*Mcosine(N,:))/Ns;
            fft_2 = fft(dfI_dqL.*Msine(N,:))/Ns;
            
            dRdXt_nl(I,:,1,L,N,1) = 2*real(fft_1(1:H+1));
            dRdXt_nl(I,:,1,L,N,2) = 2*real(fft_2(1:H+1));
            
            dRdXt_nl(I,:,2,L,N,1) = -2*imag(fft_1(1:H+1));
            dRdXt_nl(I,:,2,L,N,2) = -2*imag(fft_2(1:H+1));
        end
    end
end
dRdXt_nl(:,1,1,:,:,:) =  dRdXt_nl(:,1,1,:,:,:)/2; %divide by 2 cosine term of first harmonic

%sum linear and nonlinear contribution to gradient
dRdXt = dRdXt + dRdXt_nl;


%% dRdOmt                                                           
% derivative with respect to omega (we assume that nonlinear terms are 
% depending on displacements only)
dRdOmt = zeros(n,H+1,2);

Jmat = repmat((0:1:H),n,1); %jacobian matrix
Jmat2 = Jmat.^2;
dRdOmt(:,:,1) = -2*Om*M*squeeze(Xt(:,:,1)).*Jmat2 + D*squeeze(Xt(:,:,2)).*Jmat;
dRdOmt(:,:,2) = -2*Om*M*squeeze(Xt(:,:,2)).*Jmat2 - D*squeeze(Xt(:,:,1)).*Jmat;


%% dRdPt                                                            

dRdP_nlt = zeros(n,H+1,2,m);%derivative of nonlinear part of res
dRdP_lt = zeros(n,H+1,2,m); %derivative of linear part of res

%fill derivatives of nonlinear terms
for I = 1:n
    for N = 1:m
        fft_p = fft(squeeze(dfdp(I,N,:))')/Ns; %column vector
        
        dRdP_nlt(I,:,1,N) = 2*real(fft_p(1:H+1));
        dRdP_nlt(I,:,2,N) = -2*imag(fft_p(1:H+1));
    end
end
dRdP_nlt(:,1,1,:) = dRdP_nlt(:,1,1,:)/2; %divide static term by 2

%fill derivatives of linear terms
for J = 0:H
          
   dRdP_lt(:,J+1,1,:) =  - J^2*Om^2*double(ttv(dMdp,Xt(:,J+1,1),2)) +...
                       + J*Om*double(ttv(dDdp,Xt(:,J+1,2),2))  +...
                       + double(ttv(dKdp,Xt(:,J+1,1),2));
   
   dRdP_lt(:,J+1,2,:) =  - J^2*Om^2*double(ttv(dMdp,Xt(:,J+1,2),2)) +...
                       - J*Om*double(ttv(dDdp,Xt(:,J+1,1),2))  +...
                       + double(ttv(dKdp,Xt(:,J+1,1),2));  
end
    
%sum linear and nonlinear contributions to total derivative
dRdPt = dRdP_nlt + dRdP_lt; 

%% dRdX2t                                                           
%second order derivative of the linear part of the residual is null.
%consider only nonlinear part

dRdX2t = zeros(n,H+1,2,n,H+1,2,n,H+1,2);

% dfdq2 = fun_dfdqdq(q); %last dimension reserved for time sample

for I = 1:n
    for L = 1:n
        for O = 1:n
            
            dfI_dq2LO = squeeze(dfdq2(I,L,O,:))'; %vector of time samples of the derivative
            
            for N = 1:H+1
                for P = 1:H+1
                    
                    fft_11 = fft(dfI_dq2LO.*Mcosine(P,:).*Mcosine(N,:))/Ns;
                    fft_12 = fft(dfI_dq2LO.*Mcosine(P,:).*Msine(N,:))/Ns;
                    fft_21 = fft(dfI_dq2LO.*Msine(P,:).*Mcosine(N,:))/Ns;
                    fft_22 = fft(dfI_dq2LO.*Msine(P,:).*Msine(N,:))/Ns;
                    
                    dRdX2t(I,:,1,L,N,1,O,P,1) = 2*real(fft_11(1:H+1));
                    dRdX2t(I,:,1,L,N,2,O,P,1) = 2*real(fft_12(1:H+1));
                    dRdX2t(I,:,1,L,N,1,O,P,2) = 2*real(fft_21(1:H+1));
                    dRdX2t(I,:,1,L,N,2,O,P,2) = 2*real(fft_22(1:H+1));
                    
                    dRdX2t(I,:,2,L,N,1,O,P,1) = -2*imag(fft_11(1:H+1));
                    dRdX2t(I,:,2,L,N,2,O,P,1) = -2*imag(fft_12(1:H+1));
                    dRdX2t(I,:,2,L,N,1,O,P,2) = -2*imag(fft_21(1:H+1));
                    dRdX2t(I,:,2,L,N,2,O,P,2) = -2*imag(fft_22(1:H+1));
                    
                end
            end
            
        end
    end
end
dRdX2t(:,1,1,:,:,:,:,:,:) = dRdX2t(:,1,1,:,:,:,:,:,:)/2;


%% dRdXdOmt, dRdOm2t                                                
% (assumption that nonlinear forces are function only of displacements) ->
% -> residual of nonlinear function is independent from omega
dRdXdOmt = zeros(n,H+1,2,n,H+1,2);

for J = 1:H+1
    
    %derive cosine terms of ft of residual
    dRdXdOmt(:,J,1,:,J,1) = -2*M*Om*(J-1)^2; %derive w.r.t. A
    dRdXdOmt(:,J,1,:,J,2) = D*(J-1); %derive w.r.t. B
    
    %derive sine terms of ft of residual
    dRdXdOmt(:,J,2,:,J,1) = -D*(J-1);
    dRdXdOmt(:,J,2,:,J,2) = -2*M*Om*(J-1)^2;

end

dRdOm2t = zeros(n,H+1,2);

dRdOm2t(:,:,1) = -2*M*squeeze(Xt(:,:,1)).*Jmat2;
dRdOm2t(:,:,2) = -2*M*squeeze(Xt(:,:,2)).*Jmat2;


%% dRdXdPt, dRdPdOmt                                                

dRdXdP_nlt = zeros(n,H+1,2,n,H+1,2,m); %derivative of linear contribution of residual
dRdXdP_lt = zeros(n,H+1,2,n,H+1,2,m); %derivative of nonlin contribution of residual

for I = 1:n
    
    for L = 1:n
        
        for O = 1:m
            
        dfI_dqLdpO = squeeze(dfdqdp(I,L,O,:))'; %vector of function samples in time domain, last dimension is time sample
        
            for N = 1:H+1 %for every harmonic except the first one

                fft_1 = fft(dfI_dqLdpO.*Mcosine(N,:))/Ns;
                fft_2 = fft(dfI_dqLdpO.*Msine(N,:))/Ns;

                dRdXdP_nlt(I,:,1,L,N,1,O) = 2*real(fft_1(1:H+1));
                dRdXdP_nlt(I,:,1,L,N,2,O) = 2*real(fft_2(1:H+1));

                dRdXdP_nlt(I,:,2,L,N,1,O) = -2*imag(fft_1(1:H+1));
                dRdXdP_nlt(I,:,2,L,N,2,O) = -2*imag(fft_2(1:H+1));
            end
            
        end
    end
end
dRdXdP_nlt(:,1,1,:,:,:) =  dRdXdP_nlt(:,1,1,:,:,:)/2; %divide by 2 cosine term of first harmonic

%derivative of linear part of residual
dMdp = double(dMdp);
dKdp = double(dKdp);
dDdp = double(dDdp);

for J = 0:H
    
  dRdXdP_lt(:,J+1,1,:,J+1,1,:) = - dMdp*Om^2*J^2 + dKdp;
  dRdXdP_lt(:,J+1,1,:,J+1,2,:) = + dDdp*Om*J;
  
  dRdXdP_lt(:,J+1,2,:,J+1,1,:) = - dDdp*Om*J;
  dRdXdP_lt(:,J+1,2,:,J+1,2,:) = - dMdp*Om^2*J^2 + dKdp;
      
end

dRdXdPt = dRdXdP_nlt + dRdXdP_lt; %sum linear and nonlinear contribution to derivatives

%dRdOm
dRdPdOmt = zeros(n,H+1,2,m); 

dMdp = tensor(dMdp);
dDdp = tensor(dDdp);

for J = 0:H
    
    dRdPdOmt(:,J+1,1,:) = -2*Om*J^2*double(ttv(dMdp,Xt(:,J+1,1),2)) +...
                          + J*double(ttv(dDdp,Xt(:,J+1,2),2));
                      
    dRdPdOmt(:,J+1,2,:) = -2*Om*J^2*double(ttv(dMdp,Xt(:,J+1,2),2)) +...
                          - J*double(ttv(dDdp,Xt(:,J+1,1),2));    
end

% %we assume that the term from nonlinear term of residual is null. 
% this is true if nonlin parameter dependent terms depend on displacements 
% (and not velocities). If non linear forces are depending from speeds or
% acceleration this does not hold. 


%% dRdP2t                                                           

dRdP2t = zeros(n,H+1,2,m,m);

%fill dRdP2t tensor
for I = 1:n
    for N = 1:m
        for M = 1:m
            fft_p2 = fft(squeeze(dfdp2(I,N,M,:))')/Ns; %column vector

            dRdP2t(I,:,1,N,M) = 2*real(fft_p2(1:H+1));
            dRdP2t(I,:,2,N,M) = -2*imag(fft_p2(1:H+1));
        end
    end
end
dRdP2t(:,1,1,:,:) = dRdP2t(:,1,1,:,:)/2; %divide static term by 2


%% Convert to Nlvib representation                                  

% index of cosine - sine terms saved in matrices
ICm = zeros(n,H+1); %index cosine
I0 = (1:1:n)'; %index cosine 0roth harm

ICm(:,1) = I0;

for h = 2:H+1
    ICm(:,h) = ICm(:,h-1) + 2*n;
end

ISm = ICm + n; %index sine, rows--> dof, column--> harm number

% dRdX_____________________________________________________________________
dRdX = zeros(2*n*(H+1), n*2*(H+1)+1); %with extra dimension
for h1 = 1:H+1
    for h2 = 1:H+1
        dRdX( ICm(:,h1),  ICm(:,h2) ) = dRdXt(:,h1,1,:,h2,1);
        dRdX( ICm(:,h1),  ISm(:,h2) ) = dRdXt(:,h1,1,:,h2,2);
        dRdX( ISm(:,h1),  ISm(:,h2) ) = dRdXt(:,h1,2,:,h2,2);
        dRdX( ISm(:,h1),  ICm(:,h2) ) = dRdXt(:,h1,2,:,h2,1);  
    end    
    dRdX(ICm(:,h1),end) = dRdOmt(:,h1,1);
    dRdX(ISm(:,h1),end) = dRdOmt(:,h1,2); 
end

%eliminate derivative (that has no sense) w.r.t. HB coefficient of sine, harmonic 0
dRdX(:, n+1:2*n) = [];
dRdX(n+1:2*n, :) = [];

% dRdP_____________________________________________________________________
dRdP = zeros(2*n*(H+1),m); %with extra dimension

for h1 = 1:H+1
    dRdP( ICm(:,h1), : ) = dRdPt(:,h1,1,:);
    dRdP( ISm(:,h1), : ) = dRdPt(:,h1,2,:);   
end
dRdP(n+1:2*n,:) = [];

% dRdX2____________________________________________________________________
dRdX2 = zeros(2*n*(H+1), 2*n*(H+1) +1, 2*n*(H+1) +1); %with extra dimension
for h1 = 1:H+1
    for h2 = 1:H+1
        for h3 = 1:H+1
            dRdX2( ICm(:,h1), ICm(:,h2), ICm(:,h3) ) = dRdX2t(:,h1,1,:,h2,1,:,h3,1);
            dRdX2( ICm(:,h1), ICm(:,h2), ISm(:,h3) ) = dRdX2t(:,h1,1,:,h2,1,:,h3,2);
            
            dRdX2( ICm(:,h1), ISm(:,h2), ICm(:,h3) ) = dRdX2t(:,h1,1,:,h2,2,:,h3,1);
            dRdX2( ICm(:,h1), ISm(:,h2), ISm(:,h3) ) = dRdX2t(:,h1,1,:,h2,2,:,h3,2);
            
            dRdX2( ISm(:,h1), ISm(:,h2), ICm(:,h3) ) = dRdX2t(:,h1,2,:,h2,2,:,h3,1);
            dRdX2( ISm(:,h1), ISm(:,h2), ISm(:,h3) ) = dRdX2t(:,h1,2,:,h2,2,:,h3,2);
            
            dRdX2( ISm(:,h1), ICm(:,h2), ICm(:,h3) ) = dRdX2t(:,h1,2,:,h2,1,:,h3,1);
            dRdX2( ISm(:,h1), ICm(:,h2), ISm(:,h3) ) = dRdX2t(:,h1,2,:,h2,1,:,h3,2);
        end
        
        dRdX2( ICm(:,h1), ICm(:,h2), end ) = dRdXdOmt(:,h1,1,:,h2,1);
        dRdX2( ICm(:,h1), ISm(:,h2), end ) = dRdXdOmt(:,h1,1,:,h2,2);
        dRdX2( ISm(:,h1), ISm(:,h2), end ) = dRdXdOmt(:,h1,2,:,h2,2);
        dRdX2( ISm(:,h1), ICm(:,h2), end ) = dRdXdOmt(:,h1,2,:,h2,1); 
        
        %mixed derivatives are the same for inverted derivation order
        dRdX2( ICm(:,h1), end, ICm(:,h2) ) = dRdXdOmt(:,h1,1,:,h2,1);
        dRdX2( ICm(:,h1), end, ISm(:,h2) ) = dRdXdOmt(:,h1,1,:,h2,2);
        dRdX2( ISm(:,h1), end, ISm(:,h2) ) = dRdXdOmt(:,h1,2,:,h2,2);
        dRdX2( ISm(:,h1), end, ICm(:,h2) ) = dRdXdOmt(:,h1,2,:,h2,1); 
    end
    
    dRdX2( ICm(:,h1), end, end) = dRdOm2t(:,h1,1);
    dRdX2( ISm(:,h1), end, end) = dRdOm2t(:,h1,2);
end

%remove derivatives with respect to sine coeff of harmonic 0.
dRdX2(:, n+1:2*n, :) = [];
dRdX2(n+1:2*n, :, :) = [];
dRdX2(:, :, n+1:2*n) = [];

% dRdXdP___________________________________________________________________
dRdXdP = zeros(n*2*(H+1),n*2*(H+1)+1,m);
for h1 = 1:H+1
    for h2 = 1:H+1
        dRdXdP( ICm(:,h1),  ICm(:,h2), : ) = dRdXdPt(:,h1,1,:,h2,1,:);
        dRdXdP( ICm(:,h1),  ISm(:,h2), : ) = dRdXdPt(:,h1,1,:,h2,2,:);
        dRdXdP( ISm(:,h1),  ISm(:,h2), : ) = dRdXdPt(:,h1,2,:,h2,2,:);
        dRdXdP( ISm(:,h1),  ICm(:,h2), : ) = dRdXdPt(:,h1,2,:,h2,1,:);  
    end    
    dRdXdP(ICm(:,h1),end,:) = dRdPdOmt(:,h1,1,:);
    dRdXdP(ISm(:,h1),end,:) = dRdPdOmt(:,h1,2,:); 
end

%eliminate derivative (that has no sense) w.r.t. HB coefficient of sine, harmonic 0
dRdXdP(:, n+1:2*n,:) = [];
dRdXdP(n+1:2*n, :,:) = [];

% dRdP2 ___________________________________________________________________
dRdP2 = zeros(2*n*(H+1),m,m); %with extra dimension
for h1 = 1:H+1
    dRdP2( ICm(:,h1), :, :) = dRdP2t(:,h1,1,:,:);
    dRdP2( ISm(:,h1), :, :) = dRdP2t(:,h1,2,:,:);   
end
dRdP2(n+1:2*n,:,:) = [];


%% Convert to struct array                                          

HBderivatives.dRdX = dRdX;
HBderivatives.dRdX2 = dRdX2;
HBderivatives.dRdP = dRdP;
HBderivatives.dRdP2 = dRdP2;
HBderivatives.dRdXdP = dRdXdP;


end


%% Auxiliary function                                               
function out = extend_der_in_time(q,functionHandle)

% this function is created to have stored in a single arrays all the 
% derivatives evaluated for different time snapshots of q (displacements).
% 
% Inputs: 
%           (1) functionHandle: containts the function handle to evaluate
%                               internal forces derivatives as a function
%                               of q_ii, that is the vector of
%                               generalized displacements at time instant
%                               ii. Output of fucntionHandle must contain
%                               the following fieds: dfdp, dfdq, dfdq2,
%                               dfdp2, dfdqdp.
%
%           (2) q:              nxn_s matrix containing time snapshots of
%                               generalized displacements. q(:,ii) are
%                               displacements at time ii.
% Outputs: 
%           (1) out:            for each field of the output of
%                               @(q)functionHandle an extra dimension 
%                               corresponding to time instant is added.
%                               eg. dfdp(:,:,ii) = funcHandle(q(:,ii).dfdp)


n_s = size(q,2); %number of samples in time domain

q0 = q(:,1);

der  = functionHandle(q0);

% check output
field_names = {'dfdp','dfdq','dfdq2','dfdp2','dfdqdp'};
check_input = isfield(der,field_names);

if ~isempty(find(check_input == 0))
 error('output of @(q) funHandle(q) must contain the following fields: dfdp, dfdq, dfdq2, dfdp2, dfdqdp');
end

n = size(der.dfdp,1);
m = size(der.dfdp,2);

dfdp = zeros(n,m,n_s);
dfdq = zeros(n,n,n_s);
dfdq2 = zeros(n,n,n,n_s);
dfdp2 = zeros(n,m,m,n_s);
dfdqdp = zeros(n,n,m,n_s);

for ii = 1:n_s
    
    q_ii = q(:,ii);
    
    der = functionHandle(q_ii);
    
    dfdp(:,:,ii) = der.dfdp;
    dfdq(:,:,ii) = der.dfdq;
    dfdq2(:,:,:,ii) = der.dfdq2;
    dfdp2(:,:,:,ii) = der.dfdp2;
    dfdqdp(:,:,:,ii) = der.dfdqdp;
    
end
    
out.dfdp = dfdp;
out.dfdq = dfdq;
out.dfdq2 = dfdq2;
out.dfdp2 = dfdp2;
out.dfdqdp = dfdqdp;

end


