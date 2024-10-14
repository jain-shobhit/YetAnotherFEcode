%-------------------------------------------------------------
%
%    Copyright (C) 2007 Krister Svanberg
%
%    This file, gctoymain.m, is part of GCMMA-MMA-code.
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
%  Version July 2007.
%
%  This file contains a main program for using GCMMA to solve
%  a problem defined by the users files gctoyinit.m
%  (which must be run before gctoymain.m), toy1.m and toy2.m.
%
%%%% If outeriter=0, the user should now calculate function values
%%%% and gradients of the objective- and constraint functions at xval.
%%%% The results should be put in f0val, df0dx, fval and dfdx:
if outeriter < 0.5
  [f0val,df0dx,fval,dfdx] = toy2(xval);
  innerit=0;
  outvector1 = [outeriter innerit xval']
  outvector2 = [f0val fval']
end
%
%%%% The outer iterations start:
kktnorm = kkttol+10;
outit = 0;
while kktnorm > kkttol & outit < maxoutit
  outit   = outit+1;
  outeriter = outeriter+1;
%%%% The parameters low, upp, raa0 and raa are calculated:
  [low,upp,raa0,raa] = ...
  asymp(outeriter,n,xval,xold1,xold2,xmin,xmax,low,upp, ...
        raa0,raa,raa0eps,raaeps,df0dx,dfdx);
%%%% The GCMMA subproblem is solved at the point xval:
  [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp] = ...
  gcmmasub(m,n,outeriter,epsimin,xval,xmin,xmax,low,upp, ...
           raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d);
%%%% The user should now calculate function values (no gradients)
%%%% of the objective- and constraint functions at the point xmma
%%%% ( = the optimal solution of the subproblem).
%%%% The results should be put in f0valnew and fvalnew.
  [f0valnew,fvalnew] = toy1(xmma);
%%%% It is checked if the approximations are conservative:
  [conserv] = concheck(m,epsimin,f0app,f0valnew,fapp,fvalnew);
%%%% While the approximations are non-conservative (conserv=0),
%%%% repeated inner iterations are made:
  innerit=0;
  if conserv == 0
    while conserv == 0 & innerit <= 15
      innerit = innerit+1;
%%%% New values on the parameters raa0 and raa are calculated:
      [raa0,raa] = ...
      raaupdate(xmma,xval,xmin,xmax,low,upp,f0valnew,fvalnew, ...
                f0app,fapp,raa0,raa,raa0eps,raaeps,epsimin);
%%%% The GCMMA subproblem is solved with these new raa0 and raa:
      [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp] = ...
      gcmmasub(m,n,outeriter,epsimin,xval,xmin,xmax,low,upp, ...
               raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d);
%%%% The user should now calculate function values (no gradients)
%%%% of the objective- and constraint functions at the point xmma
%%%% ( = the optimal solution of the subproblem).
%%%% The results should be put in f0valnew and fvalnew:
      [f0valnew,fvalnew] = toy1(xmma);
%%%% It is checked if the approximations have become conservative:
      [conserv] = concheck(m,epsimin,f0app,f0valnew,fapp,fvalnew);
    end
  end
%%%% No more inner iterations. Some vectors are updated:
  xold2 = xold1;
  xold1 = xval;
  xval  = xmma;
%%%% The user should now calculate function values and gradients
%%%% of the objective- and constraint functions at xval.
%%%% The results should be put in f0val, df0dx, fval and dfdx:
  [f0val,df0dx,fval,dfdx] = toy2(xval);
%%%% The residual vector of the KKT conditions is calculated:
  [residu,kktnorm,residumax] = ...
  kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s, ...
           xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);
  outvector1 = [outeriter innerit xval']
  outvector2 = [f0val fval']
end
%---------------------------------------------------------------------

