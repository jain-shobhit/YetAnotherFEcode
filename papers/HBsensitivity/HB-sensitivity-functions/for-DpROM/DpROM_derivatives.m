% DESCRIPTION:
% DpROM_derivatives provides as output the derivatives with respect to
% displacements and parameters of the internal forces of DpROM model 
% described in ref*.
% The e.o.m. for this FE model is:
% Mq_dd + Dq_d + K(P)q + f_nl(q,p) = f_ex, in which p is the parameter
% vector of defects amplitudes. For more detail see [1].
%
% INPUTS: 
% (1) q:                vector of time domain displacements
% (2) tensors_DpROM:    structure array containing reduced tensors of Fe
%                       model described in [1]
% (3) varargin{1}:      Set to 1 if the field of defect displacements
%                       resulst in an isochoric transformation. This allows
%                       to reduce computational effort. If this input is
%                       not inserted the funciton computes derivatives as
%                       if the defect is not isochoric.
%
% OUTPUTS:
% (1) der:              strucure array containing first order, second order
%                       and mixed derivatives with respect to displacement
%                       vector q and parameter vector p. These derivatives
%                       are stored in the following fields:
%                       'dfdp','dfdq',dfdq2','dfdp2','dfdqdp','dMdp'
% 
% REFERENCES:
% [1] Marconi,Tiso,Quadrelli,Braghin 'A higher-order parametric nonlinear
% reduced-order model for imperfect structures using Neumann expansion',
% Nonlinear Dynamics 104(4), http://doi.org/10.1007/s11071-021-06496-y
%
% Author: Alexander Saccani, Msc in mechanical engineering
% University: Politecnico di Milano
% Created:  09/2021

function der = DpROM_derivatives(q, tensors_DpROM, varargin)

% Set input
if nargin < 3
    isIsochoric = 0;    % default option computes derivatives considering a 
                        % volume change (stemming from defect) 
else
    isIsochoric = varargin{1};  % if nargin is set to 1, volume is 
                                % considered constant when applying defect
end


% *************************************************************************
% populate tensors
% *************************************************************************

Q2n = tensors_DpROM.Q2n{1,1}; 
Q3n = tensors_DpROM.Q3n{1,1};
Q4n = tensors_DpROM.Q4n{1,1};
Q3d = tensors_DpROM.Q3d{1,1};
Q4d = tensors_DpROM.Q4d{1,1};
Q5d = tensors_DpROM.Q5d{1,1};
Q4dd = tensors_DpROM.Q4dd{1,1};
Q5dd = tensors_DpROM.Q5dd{1,1};
Q6dd = tensors_DpROM.Q6dd{1,1};

M = tensors_DpROM.M; %mass matrices, 

% length of parameter vector
if length(size(Q3d)) == 2 
    m = 1;
else
    m = size(Q3d,3);
end

%create tensors useful for non isochoric transformation
if isIsochoric ~= 1
    %initialize tensors
    Q2n_ext = tensor(zeros([size(Q2n),m])); %'ext' stands for 'extended'
    Q3n_ext = tensor(zeros([size(Q3n),m]));
    Q4n_ext = tensor(zeros([size(Q4n),m]));
    Q3d_ext = tensor(zeros([size(Q3d),m]));
    Q4d_ext = tensor(zeros([size(Q4d),m]));
    Q5d_ext = tensor(zeros([size(Q5d),m]));
    %fill useful tensors for non isochoric transformation
    for ii = 1:m      
        Q2n_ext(:,:,ii)         = tensors_DpROM.Q2n{1, ii+1};
        Q3n_ext(:,:,:,ii)       = tensors_DpROM.Q3n{1, ii+1};
        Q4n_ext(:,:,:,:,ii)     = tensors_DpROM.Q4n{1, ii+1};
        Q3d_ext(:,:,:,ii)       = tensors_DpROM.Q3d{1, ii+1};  
        Q4d_ext(:,:,:,:,ii)     = tensors_DpROM.Q4d{1, ii+1};
        Q5d_ext(:,:,:,:,:,ii)   = tensors_DpROM.Q5d{1, ii+1};
    end 
    %adjust tensor if m == 1
    if m == 1
        Q2n_ext = squeeze(Q2n_ext);
        Q3n_ext = squeeze(Q3n_ext);
        Q4n_ext = squeeze(Q4n_ext);
        Q3d_ext = squeeze(Q3d_ext);
        Q4d_ext = squeeze(Q4d_ext);
        Q5d_ext = squeeze(Q5d_ext);
    end
end


% *************************************************************************
% Compute derivatives of internal forces and mass matrix
% *************************************************************************

if m==1
    dfdq =  ttv(Q3n,q,3) + ttv(Q3n,q,2) + ttv(Q4n,{q,q},[3,4]) +...
        ttv(Q4n,{q,q},[2,3])+ttv(Q4n,{q,q},[2,4]); %without linear part

    dfdp =  ttv(Q3d,q,2) + ttv(Q4d,{q,q},[2 3]) + ...
        +ttv(Q5d,{q,q,q},[2 3 4]); % with also the linear part

    dfdq2 = Q3n + permute(Q3n,[1,3,2]) + ttv(Q4n,q,4) +...
        + ttv(Q4n,q,3) + permute(ttv(Q4n,q,3)+ttv(Q4n,q,2),[1,3,2]) +...
        permute(ttv(Q4n,q,4),[1,3,2])+ttv(Q4n,q,2); % without linear part

    dfdp2 = 2*ttv(Q4dd,q,2) + 2*ttv(Q5dd,{q,q},[2,3]) + ...
        2*ttv(Q6dd,{q,q,q},[2,3,4]);

    dfdpdq = Q3d + ttv(Q4d,q,3) + ttv(Q4d,q,2) + ...
        ttv(Q5d,{q,q},[3,4]) + ttv(Q5d,{q,q},[2,4]) + ...
        ttv(Q5d,{q,q},[2,3]);

    dfdqdp = dfdpdq;

    if isIsochoric ~= 1 % additional terms for not-isochoric defects
        dfdp =  dfdp + ttv(Q2n_ext,q,2) + ttv(Q3n_ext,{q,q},[2,3]) +...
            + ttv(Q4n_ext,{q,q,q},[2,3,4]);

        dfdp2 = dfdp2 + ttv(Q3d_ext,q,2) + ttv(Q3d_ext,q,2)+ ...
            + ttv( 2*Q4d_ext, {q,q}, [2,3]) + ...
            + ttv( 2*Q5d_ext, {q,q,q}, [2,3,4]);

        dfdpdq = dfdpdq + Q2n_ext + ...
            ttv(Q3n_ext,q,3) + ttv(Q3n_ext,q,2) + ...
            ttv(Q4n_ext,{q,q},[3,4])+ttv(Q4n_ext,{q,q},[2,4]) + ...
            ttv(Q4n_ext,{q,q},[2,3]);

        dfdqdp = dfdpdq; %puoi togliere se modifichi codice
    end
else 
    dfdq =  ttv(Q3n,q,3) + ttv(Q3n,q,2) + ttv(Q4n,{q,q},[3,4]) +...
        ttv(Q4n,{q,q},[2,3]) + ttv(Q4n,{q,q},[2,4]); % without linear part

    dfdp =  ttv(Q3d,q,2) + ttv(Q4d,{q,q},[2 3]) + ...
        +ttv(Q5d,{q,q,q},[2 3 4]); % with also the linear part

    dfdq2 = Q3n + permute(Q3n,[1,3,2]) + ttv(Q4n,q,4) +...
        + ttv(Q4n,q,3) + permute(ttv(Q4n,q,3)+ttv(Q4n,q,2),[1,3,2]) +...
        + permute(ttv(Q4n,q,4),[1,3,2])+ttv(Q4n,q,2); % without linear part

    dfdp2 = permute(ttv(Q4dd,q,2),[1,3,2]) + ttv(Q4dd,q,2) +...
        ttv(Q5dd,{q,q},[2,3]) + permute(ttv(Q5dd,{q,q},[2,3]),[1,3,2]) +...
        permute(ttv(Q6dd,{q,q,q},[2,3,4]),[1,3,2]) + ttv(Q6dd,{q,q,q},[2,3,4]);

    dfdpdq = permute(Q3d ,[1,3,2]) + permute(ttv(Q4d,q,3) + ...
        ttv(Q4d,q,2),[1,3,2])+ permute(ttv(Q5d,{q,q},[3,4])+...
        + ttv(Q5d,{q,q},[2,4])+ttv(Q5d,{q,q},[2,3]),[1,3,2]);

    dfdqdp = permute(dfdpdq,[1,3,2]); %puoi togliere se modifichi codice

    if isIsochoric ~= 1 % additional terms for not-isochoric defects
        dfdp =  dfdp + ttv(Q2n_ext,q,2) + ttv(Q3n_ext,{q,q},[2,3]) +...
            + ttv(Q4n_ext,{q,q,q},[2,3,4]); 

        dfdp2 = dfdp2 + ttv(Q3d_ext,q,2) + permute(ttv(Q3d_ext,q,2),[1,3,2])+...
            + ttv( Q4d_ext + permute(Q4d_ext,[1,2,3,5,4]), {q,q}, [2,3]) + ...
            + ttv( Q5d_ext + permute(Q5d_ext,[1,2,3,4,6,5]), {q,q,q}, [2,3,4]);

        dfdpdq = dfdpdq + permute(Q2n_ext,[1,3,2]) + ...
            permute(ttv(Q3n_ext,q,3) + ttv(Q3n_ext,q,2),[1,3,2]) + ...
            permute( ttv(Q4n_ext,{q,q},[3,4])+ttv(Q4n_ext,{q,q},[2,4]) + ...
            ttv(Q4n_ext,{q,q},[2,3]), [1,3,2]);

        dfdqdp = permute(dfdpdq,[1,3,2]); %puoi togliere se modifichi codice
    end
    
    if isIsochoric ~= 1
        
        n = size(M{1},1);
        dMdp = zeros(n,n,m);
        
        for J = 1:m
            dMdp(:,:,J) = M{J+1};
        end
        
    end 
end 


% *************************************************************************
% Generate output
% *************************************************************************
der.dfdp = double(dfdp);
der.dfdq = double(dfdq);
der.dfdq2 = double(dfdq2);
der.dfdp2 = double(dfdp2);
der.dfdqdp = double(dfdqdp);

if isIsochoric ~= 1
    der.dMdp = dMdp;
end

end