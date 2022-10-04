% nlvib_StabilityAnalysis
%
% [stable, mucrit, X_shoot, a, a_rms] = ...
%       nlvib_StabilityAnalysis(mySystem, X, method, H, Nt, Np)
%
% Description: compute stability for HB and shooting frequency response
% curves obtained through NLvib.
%
% INPUTS
%   - mySystem: system used in NLvib
%   - X: frequency response solution from NLvib
%   - method: shooting/HB (choose the one used to compute X)
%   - H: number of harmonics used for the HB (if method is HB), and/or 
%     number harmonics considered in a_rms.
%   - Nt*: number of time samples per period (*optional, default is 2^8)
%   - Np*: seek for Np-period solutions (*optional, default is 1)
% OUTPUTS
%   - stable: vector of 0 (unstable) and 1 (stable)
%   - mucrit: vector containing the largest magnitude monodromy matrix
%     eigenvalue
%   - X_shoot: returns X in NLvib's shooting format (if method is already
%     shooting, then X_shoot=X)
%   - a: matrix containing, by columns, the first hamonic magnitude for each
%     dof of the system
%   - a_rms: same as "a", but with the RMS value computed for the first H
%     harmonics
%
% Author: Jacopo Marconi, PhD, Politecnico di Milano
% Last modified: 5/8/2022

function [stable, mucrit, X_shoot, a, a_rms] = nlvib_StabilityAnalysis(mySystem, X, method, H, Nt, Np)

if nargin < 5
    Nt = 2^8;  % Number of time samples per period
    Np = 1;    % we seek period-one solutions
elseif nargin < 6
    Np = 1;    % we seek period-one solutions
end

if strcmpi(method,'shooting')
    X_shoot = X;
elseif strcmpi(method,'HB')
    % convert HB solutions in shooting solutions
    X_HB = X;
    Q  = X_HB(1:end-1,:);
    Om = X_HB(end,:);
    n = (size(X_HB,1) - 1) / (2*H + 1); % number of dofs (assuming H harmonics)
    % RECONSTRUCT SOLUTION AT TIME=0s
    I0 = 1:n; % 0th order term indexes
    IC = n+repmat(1:n,1,H)+n*kron(0:2:2*(H-1),ones(1,n)); % cos indexes
    IS = IC+n; % sine indexes
    kk = kron((1:H)', ones(n,1));
    q = zeros(n, length(Om));
    u = zeros(n, length(Om));
    for ww = 1 : length(Om)
        Omega = Om(ww);
        t = 2*pi/Omega;
        % add harmonics
        Q0  = Q(I0, ww);
        Qre = Q(IC, ww);
        Qim = Q(IS, ww);
        
        % EXPONENTIAL representation
%         C = (Qre - 1i*Qim).*exp(1i*kk*Omega*t); % watch out for the minus sign!
%         q(:,ww) = real(sum(reshape([Q0; C], n, H+1), 2));
%         u(:,ww) = real(sum(reshape((1i*kk*Omega).*C, n, H), 2)) / Omega;        
        % SINE-COSINE representation
        CC = Qre.*cos(kk*Omega*t) + Qim.*sin(kk*Omega*t);
        CS =-Qre.*sin(kk*Omega*t) + Qim.*cos(kk*Omega*t);
        q(:,ww) = real(sum(reshape([Q0; CC], n, H+1), 2));
        u(:,ww) = real(sum(reshape((kk*Omega).*CS, n, H),2)) / Omega; % /Omega is the normalization used in X_shoot
    end
    X_shoot = [q; u; Om];
end

% Interpret solver output
Om = X_shoot(end,:);
Ys = X_shoot(1:end-1,:);

% Calculate amplitude also for the results of the shooting method, and
% determine the asymptotic stability according to Floquet theory
h = 1;
a = zeros(length(Om), size(Ys, 1));
a_rms = a;
stable  = zeros(length(Om), 1);
mucrit  = zeros(length(Om), 1);
tic;
for ii = 1:length(Om)
    % Evaluate solution and monodromy matrix
    [~, ~, ~, Y, dye_dys] = shooting_residual(X_shoot(:,ii),...
        mySystem, Nt, Np, 'FRF');
    
    % Determine fundamental harmonic magnitude
    Qc = fft(Y)/Nt;
    a(ii,:) = 2*abs( Qc(1+h, :) );
    a_rms(ii, :) = sqrt(sum(abs(Qc(1:H+1, :)).^2))/sqrt(2)*2;
    
    % Determine stability in accordance with Floquet theory: a periodic
    % solution is stable, if all eigenvalues of the monodromy matrix remain
    % within the unit circle in the complex plane
    mucrit(ii) = eigs(dye_dys, 1, 'lm'); % leading Floquet multiplier
    stable(ii) = abs(mucrit(ii))<=1;      % allow for some tolerance
end
disp(['A posteriori Floquet stability analysis required ' ...
    num2str(toc) ' s.']);

% figure;
% ii = 1;
% a = a_shoot;
% as = a(:,ii);
% as(stable==0) = NaN;
% au = a(:,ii);
% au(stable==1) = NaN;
% plot(Om, a(:,ii), 'g')
% hold on
% plot(Om, as, 'b.')
% plot(Om, au, 'r.')
% title(['dof #' num2str(ii) ' (' method ')'])
% xlabel('\Omega','FontSize',20)
% ylabel(['|Q_' num2str(h) '|'],'FontSize',20)
% grid on
% drawnow

