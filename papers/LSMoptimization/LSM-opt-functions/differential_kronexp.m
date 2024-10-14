
function dpx = differential_kronexp(p,ii)
% example:
% p = [p1 p2]'
% ii = 3
% px = kronexp(p,ii)
%
%        [   p1^3]       [ 3*p1^2,       0]
%        [p1^2*p2]       [2*p1*p2,    p1^2]
%        [p1^2*p2]       [2*p1*p2,    p1^2]
%        [p1*p2^2]       [   p2^2, 2*p1*p2]
%   px = [p1^2*p2] ----> [2*p1*p2,    p1^2] = dpx
%        [p1*p2^2]       [   p2^2, 2*p1*p2]
%        [p1*p2^2]       [   p2^2, 2*p1*p2]
%        [   p2^3]       [      0,  3*p2^2]

M = length(p);
px = kronexp(p,ii);

% % compute derivative coefficients using symbolic toolbox __________________
% filename = sprintf('coefficients_M%d_i%d.mat',M,ii);
% if exist(filename,"file")>0
%     load(filename,'diffkronexp_coeffs')
% else
%     x = sym('p', [M,1]);
%     xx = kronexp(x,ii);
%     diffkronexp_coeffs = zeros(M^ii, M);
%     for c = 1 : M
%         diffkronexp_coeffs(:,c) = polynomialDegree(xx,x(c));
%     end
%     save(filename,'diffkronexp_coeffs')
% end
% % coefficients can be computed once and for all in the LSM class for a
% % number of expansion orders, given a vector p of size M.
% % _________________________________________________________________________

diffkronexp_coeffs = zeros(M^ii, M);
for jj = 1 : M
    a = zeros(M,1);
    a(jj) = 1;
    diffkronexp_coeffs(:,jj) = kronaddexp(a,ii);
end

dpx = diffkronexp_coeffs.*repmat(px,1,M)./repmat(p.',length(px),1);

% % check:
% dp_test = [diff(xx, x(1)) diff(xx, x(2))];
% dp-dp_test
