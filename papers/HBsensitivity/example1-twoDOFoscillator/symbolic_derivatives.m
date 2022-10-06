% DESCRIPTION                                                       
% This file can be used to generate derivatives function for sensitivity
% analysis of frequency response curves. What is necessary to do is to:
% (1) define the vectors of generalized coordinates and of parameters
% (2) write as a function of the generalized coordinates and of parameters
% the potential energy of linear (V_l) and of nonlinear (V_nl) forces.
% (3) execute the code and copy from the command window dfdq, dfdp, ...
% etc.. in a new matlab funciton script.
%
% Before running make sure symbolic matlab tool is installed.
%
% Author: Alexander Saccani, msc in mechanical engineering,  
%         Politecnico di Milano
% Created: 09/2021

% configure matlab
clear
clc
close all


%% Input                                                            
% declear generalized coordinates and parameters
syms q1 q2 k11 k13 k21 k23 m1 m2 d1 d2 'real'

% generalized coordinates vector
q = [q1,q2].';

% parameter vector
p = [k11 k13 k21 k23 m1 m2 d1 d2];

% define potential energy
V_l = 1/2*k11*q1^2 + 1/2*k21*(q2-q1)^2;  %.. of linear forces
V_nl = 1/4*k13*q1^4 + 1/4*k23*(q2-q1)^4; %.. of nonlinear forces


%% computation of derivatives                                       
n = length(q);
m = length(p);

% forces from potential energy 
f_l = sym(zeros(n,1));  % lin
f_nl = sym(zeros(n,1)); % nonlin
for ii = 1:n
    f_l(ii) = diff(V_l,q(ii));
    f_nl(ii) = diff(V_nl,q(ii));
end
f = f_l + f_nl;         % all

% initialize output 
dfdq = sym(zeros(n));
dfdq_all = sym(zeros(n));
dfdp = sym(zeros(n,m));
dfdqdp = sym(zeros(n,n,m));
dfdp2 = sym(zeros(n,m,m));
dfdq2 = sym(zeros(n,n,n));

% dfdq (only nonlinear terms in f)
for ii = 1:n
    dfdq(:,ii) = diff(f_nl, q(ii));
end

% dfdp (linear + nonlinear terms in f)
for ii = 1:m
    dfdp(:,ii)  = diff(f, p(ii));
end

% dfdq2
for ii = 1:n
    for jj = 1:n
        dfdq2(:,ii,jj) = diff(dfdq(:,ii), q(jj));
    end
end

%dfdp2 (linear + nonlinear terms in f)
for ii = 1:m
    for jj = 1:m
        dfdp2(:,ii,jj) = diff(dfdp(:,ii), p(jj));
    end
end

% dfdq_all (derivatives of total force (lin+nlin) with respect to q
for ii = 1:n
    dfdq_all(:,ii) = diff(f, q(ii));
end

% dfdqdp (linear and nonlinear terms in f
for ii = 1:n
    for jj = 1:m
        dfdqdp(:,ii,jj) = diff(dfdq_all(:,ii), p(jj));
    end
end

