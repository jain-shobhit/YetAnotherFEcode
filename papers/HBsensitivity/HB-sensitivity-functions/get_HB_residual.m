function normRes = get_HB_residual(Xsol,system,H,n,N,harm,dof,AB)
% DESCRIPTION:
% get_HB_residual function returns with a column vector containing the norm
% of the residual of harmonic balance equations for system specified in
% input. The residual euclidian norm can be computed for the full residual
% vector or can be limited only to some elements, defined by harmonic 
% coefficients and dofs in the input variables.
%
% INPUTS:
% (1) Xsol:         double array containing sets of points for which you
%                   want to evaluate the residual norm. HB coefficients in
%                   Xsol must be collected in the same order of solution
%                   output from NlVib [1]
% (2) system:       system is the same input of 'HB_residual.m' of NlVib [1]
% (3) H:            number of  harmonics
% (4) n:            number of dofs of system
% (5) N:            number of samples time to compute residual with A.F.T. 
%                   algorithm
% (6) harm:         cell array containing the harmonics that we want select
%                   from the residual vector. E.g. if harm = {1,3} the norm
%                   will be extended only to residual equations
%                   corresponding to coefficients multiplying harmonic 1
%                   and 3 of the residual. Set harm to 'all' to include
%                   all harmonics in the residual norm.
% (7) dof:          cell array containing dofs that we want to select
%                   from the residual vector. Norm will be extended only to
%                   equations representing forces equilibriums of selected
%                   dofs. Set dof to 'all' to include all dof in the
%                   residual norm.
% (8) AB:           char vector to specify wether to limit residual norm 
%                   computation to real coefficients ( set AB = 'A'), only
%                   to imaginary coefficients (set AB = 'B') or 
%                   to both real and imaginary (set AB = 'AB').
%
% OUTPUTS:
% (1) normRes:      vector array containing the norm of the residual for
%                   each point in HB coefficients space specified in Xsol.
%
% REFERENCES:
% [1]: NLvib Version 1.3, Copyright (C) 2020 Malte Krack, Johann Gross
% 
% Author: Alexander Saccani, Msc in mechanical engineering
% University: Politecnico di Milano
% Created:  09/2021

n_samples = size(Xsol,2);

I0 = (1:1:n)';
IC = n+repmat(1:n,1,H)+n*kron(0:2:2*(H-1),ones(1,n)); IS = IC+n; %cosine and sine indeces
ICm = reshape(IC,[n,H]); %row ->dof, column ->harmonic number
ISm = reshape(IS,[n,H]);

if ischar(harm)
    if strcmp(harm,'all')
        harm = 0:1:H;
    else
        error('not valid input')
    end
end

if ischar(dof)
    if strcmp(dof,'all')
        dof = 1:1:n;
    else
        error('not valid input')
    end
end

flag = 0;
logical = find( harm == 0);
if ~isempty(logical)
    harm(logical) = [];
    flag = 1;
end
 
flag1 = 0;
if isempty(harm)
    harm = 1;
    flag1 = 1;
    if flag1 && ~strcmp(AB,'A')
        error('zeroth harmonic belongs to real coefficients. Set AB to A');
    end
end

index_cos = reshape(ICm(dof,harm),[],1);
if flag == 1
    index_cos = [index_cos;I0];
end

if flag1 == 1
    index_cos = I0(dof);
end

index_sin = reshape(ISm(dof,harm),[],1);
index_cos_sin = [index_cos;index_sin];

%initialize output
normRes = zeros(n_samples,1);

for jj = 1:n_samples
    
    X = Xsol(:,jj);
    res_vec = feval(@HB_residual,X,system,H,N,'frf');
    
    switch(AB)
        case('A')
            normRes(jj) = norm(res_vec(index_cos));
        case('B')
            normRes(jj) = norm(res_vec(index_sin));
        case('AB')
            normRes(jj) = norm(res_vec(index_cos_sin));
        otherwise
            error('not valid input')
    end
            
end
        
   
end

