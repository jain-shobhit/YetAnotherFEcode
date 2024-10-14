function [fval, fgrad] = fun(x, nMasses, nSprings)
% Minimize the value of the middle mass
fval = x(2);
fgrad = zeros(1, length(x));
fgrad(2) = 1;
