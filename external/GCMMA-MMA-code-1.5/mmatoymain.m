%-------------------------------------------------------------
%
%    Copyright (C) 2009 Krister Svanberg
%
%    This file, mmatoymain.m, is part of GCMMA-MMA-code.
%    
%    GCMMA-MMA-code is free software; you can redistribute it and/or
%    modify it under the terms of the GNU General Public License as 
%    published by the Free Software Foundation; either version 3 of 
%    the License, or (at your option) any later version.
%    
%    This code is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%    
%    You should have received a copy of the GNU General Public License
%    (file COPYING) along with this file.  If not, see 
%    <http://www.gnu.org/licenses/>.
%    
%    You should have received a file README along with this file,
%    containing contact information.  If not, see
%    <http://www.smoptit.se/> or e-mail mmainfo@smoptit.se or krille@math.kth.se.
%
%------
%
%  Version September 2009.
%
%  This file contains a main program for using MMA to solve
%  a problem defined by the users files mmatoyinit.m
%  (which must be run before mmatoymain.m) and toy2.m.
%
%%%% If outeriter=0, the user should now calculate function values
%%%% and gradients of the objective- and constraint functions at xval.
%%%% The results should be put in f0val, df0dx, fval and dfdx:
%
if outeriter < 0.5
  [f0val,df0dx,fval,dfdx] = toy2(xval);
  outvector1 = [outeriter xval']
  outvector2 = [f0val fval']
end
%
%%%% The iterations start:
kktnorm = kkttol+10;
outit = 0;
while kktnorm > kkttol & outit < maxoutit
  outit   = outit+1;
  outeriter = outeriter+1;
%%%% The MMA subproblem is solved at the point xval:
  [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
  mmasub(m,n,outeriter,xval,xmin,xmax,xold1,xold2, ...
  f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);
%%%% Some vectors are updated:
  xold2 = xold1;
  xold1 = xval;
  xval  = xmma;
%%%% The user should now calculate function values and gradients
%%%% of the objective- and constraint functions at xval.
%%%% The results should be put in f0val, df0dx, fval and dfdx.
  [f0val,df0dx,fval,dfdx] = toy2(xval);
%%%% The residual vector of the KKT conditions is calculated:
  [residu,kktnorm,residumax] = ...
  kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s, ...
           xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);
  outvector1 = [outeriter xval']
  outvector2 = [f0val fval']
%
end
%---------------------------------------------------------------------
