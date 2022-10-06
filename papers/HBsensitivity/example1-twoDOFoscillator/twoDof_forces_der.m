function der = twoDof_forces_der(q,p)

%set input
q1 = q(1);
q2 = q(2);

%parameter names
k11 = p(1);
k13 = p(2);
k21 = p(3);
k23 = p(4);
m1 = p(5);
m2 = p(6);
d1 = p(7);
d2 = p(8);

% define derivatives
% dfdq (only nonlinear forces in f)
dfdq = reshape([k13.*q1.^2.*3.0+k23.*(q1-q2).^2.*3.0,k23.*(q1-q2).^2.*-...
    3.0,k23.*(q1-q2).^2.*-3.0,k23.*(q1-q2).^2.*3.0],[2,2]);

% dfdp (linear + nonlinear forces in f)
dfdp = reshape([q1,0.0,q1.^3,0.0,q1-q2,-q1+q2,(q1-q2).^3,-(q1-q2).^3,...
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[2,8]);

% dfdq2 (nonlinear forces, dfdq2 for linear forces is null)
dfdq2 = reshape([k13.*q1.*6.0+k23.*(q1.*2.0-q2.*2.0).*3.0,...
    k23.*(q1.*2.0-q2.*2.0).*-3.0,k23.*(q1.*2.0-q2.*2.0).*-3.0,...
    k23.*(q1.*2.0-q2.*2.0).*3.0,k23.*(q1.*2.0-q2.*2.0).*-3.0,...
    k23.*(q1.*2.0-q2.*2.0).*3.0,k23.*(q1.*2.0-q2.*2.0).*3.0,...
    k23.*(q1.*2.0-q2.*2.0).*-3.0],[2,2,2]);


% dfdp2 (linear + nonlinear forces)
dfdp2 = zeros(2,8,8);

% dfdpdq (linear + nonlinear forces)
dfdqdp = reshape([1.0,0.0,0.0,0.0,q1.^2.*3.0,0.0,0.0,0.0,1.0,-1.0,...
    -1.0,1.0,(q1-q2).^2.*3.0,(q1-q2).^2.*-3.0,(q1-q2).^2.*-3.0,...
    (q1-q2).^2.*3.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,...
    0.0,0.0,0.0,0.0],[2,2,8]);

% derivative of mass matrix dMdp_ijk   k -> number of parameter
% M = [m1 0; 0 m2];
dMdp = zeros(2,2,8);
dMdp(1,1,5) = 1;
dMdp(2,2,6) = 1;

% derivative of mass matrix dDdp_ijk   k -> number of parameter
% D = [d1+d2 -d2; -d2 d2]
dDdp = zeros(2,2,8);
dDdp(1,1,7) = 1;
dDdp(1,1,8) = 1;
dDdp(1,2,8) = -1;
dDdp(2,1,8) = -1;
dDdp(2,2,8) = 1;

%Notice that we didn't define dKdp. Indeed we took into account for this
%term in dfdp, dfdqdp. Anyway it is still possible to define dKdp and
%not introduce derivatives of linear stiffness terms in dfdp and dfdqdp. In
%this case add field .dKdp to output

% generate output in struct array
der.dfdq = dfdq;
der.dfdp = dfdp;
der.dfdq2 = dfdq2;
der.dfdp2 = dfdp2;
der.dfdqdp = dfdqdp;

der.dMdp = dMdp;
der.dDdp = dDdp;

end