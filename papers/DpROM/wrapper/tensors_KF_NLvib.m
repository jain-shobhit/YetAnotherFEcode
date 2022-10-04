
function [Kt,fi] = tensors_KF_NLvib(Q3,Q4,Q3t,Q4t,q)
% NB: in NLvib, the nonlinear function must contain ONLY the nonlinear
% terms (i.e. WITHOUT the linear ones)
fi = ttsv(Q3,q,-1)  + ttsv(Q4,q,-1);
Kt = ttsv(Q3t, q, -2) + ttsv(Q4t,q,-2);