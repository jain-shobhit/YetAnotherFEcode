function [fval, fgrad] = fun(x)
% Minimize the product between k2 and k3

fval = x(3)*x(4);
fgrad = [0, 0, x(4), x(3)];
