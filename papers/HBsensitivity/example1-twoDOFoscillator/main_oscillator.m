% Parameter sensitivity analysis of frequency response of a 2dof    
% spring mass damper oscillator.
% See attached pdf file for details on the system
%
% Author: Alexander Saccani, Msc in mechanical engineering, PoliMi
% Created: 09/2021

close all
clear
clc

%% Set input                                                        

% PARAMETER OF THE SYSTEM__________________________________________________
%masses:
m1 = 1; 
m2 = 1;

%springs:  k_ij : i-> dof, j-> exponent of force     
k11 = 1;
k13 = 1;
k21 = 1;
k23 = 1;

%damping:
d1 = 0.2;
d2 = 0.2;

%load:
F = 1;

% DEFINITION OF PARAMETER VECTOR AND ITS VARIATION_________________________
% with reference to the attached pdf file we define the parameter vector 
% with respect to sensitivity analysis is performed. We define its nominal
% value and its variation from the nominal.

% parameter nominal value
p = [k11, k13, k21, k23, m1, m2, d1, d2]';

% percentage parameter variation
dp_perc = [-10, -10, -10, -10, +20, +20, +30 , +30]';

% actual parameter variation
dp = p.*dp_perc/100;

% ANALYSIS SETTINGS________________________________________________________

% number of harmonics
H = 5;

% method for sensitivity analysis
method = 'normal';

% expansion order for sensitivity analysis
ord_exp = 2;


%% Nominal solution and sensitivity analysis                        

% nominal solution is computed simulating the system for parameter
% vector p. Sensitivity analysis to parameters is run in parallel with the
% nominal solution

%length of parameter vector
m = length(p);

%number of dofs
n = 2;

% CREATE SYSTEM FOR NLvib__________________________________________________
%1st cubic spring
nonlinear_elements{1} = struct('type', 'cubicspring',...
    'stiffness', k13, 'force_direction', [1;0]);

%2nd cubic spring
nonlinear_elements{2} = struct('type', 'cubicspring',...
    'stiffness', k23, 'force_direction', [1;-1]);

%generate system
system = ChainOfOscillators([m1,m2],[d1,d2,0],[k11,k21,0],...
    nonlinear_elements,[F,0]);

% EIGENVALUE ANALYSIS _____________________________________________________
%stiffness matrix
K = system.K;

%mass matrix
M = system.M;

%eigenvalue analysis
[~,om] = eigs(K, M, 2, 'SM');
[om_n,~] = sort(sqrt(diag(om)));

% SET PARAMETERS FOR ANALYSIS______________________________________________
% natural frequency of interest
foi = 1; 

% set frequency range for analysis
Om_s = .2*om_n(foi);    % start frequency
Om_e = 4*om_n(foi);     % end frequency 

% number of samples for A.F.T.
N = 3*H+1;	% for residual of nonlinear forces 
Ns = 2*N;	% for derivatives of the residual (usually this value has to be
            % greater then N)

% initial guess for continuation (linear approximation)
Q1 = (-Om_s^2*system.M + 1i*Om_s*system.D + system.K)\[F,0]';
x0 = zeros((2*H+1)*length(Q1),1);
x0(length(Q1)+(1:2*length(Q1))) = [real(Q1);-imag(Q1)]; 

% define function handle for derivatives in time domain
derFunHandle = @(q)twoDof_forces_der(q,p);

% define postprocessing function for sensitivity analysis
fun_postprocess = @(X)HB_sensitivity(X, system, derFunHandle, m, method,...
                                    ord_exp, Ns);
% define path continuation step size
ds = .01;

% other settings for continuation
Solopt = optimset(optimset(@fsolve),'Display','off',...
        'Jacobian','on');
Sopt = struct('jac','on'); % analytical Jacobian provided here

% RUN ANALYSIS_____________________________________________________________

[Xnom,SolInfo,Sol] = solve_and_continue(x0,...
    @(X) HB_residual(X,system,H,N), Om_s, Om_e, ds, ...
     Sopt, fun_postprocess, Solopt); 

% decode nominal solution
Xn = nlvib_decode(Xnom, SolInfo, Sol, 'FRF', 'HB', n, H);

% Xn is a struct array with fields 'Q0','Qre','Qim','omega'. 
% Xn.Q0{ij} --> static coefficient A0,  dof i, sample j
% Xn.Qre{ijk} --> real coefficient, dof i, sample j, harmonic k
% Xn.Qim{ijk} --> imaginary coefficient, dof i, sample j, harmonic k
% Xn.omega --> vector of ang. frequencies

% compute RMS for harmonics [0,2,3,4,5]
harm = [0, 2:H];
Xn.RMS = RMS(Xn,harm);

% create struct array with sensitivities
sens = [Sol.Sens];


%% Sensitivity based solution                                       

% response for p = p + dp computed with 1st and 2nd order sensitivities
X12 = update_response(Xnom, sens, dp);

% decode sensitivity based response
X1 = nlvib_decode(X12.lin,  0, 0, 'FRF', 'HB', n, H); % linear approx.
X2 = nlvib_decode(X12.quad, 0, 0, 'FRF', 'HB', n, H); % quadratic approx.

%compute RMS for harmonics [0,2,3,4,5]
harm = [0, 2:H];
X1.RMS = RMS(X1, harm);
X2.RMS = RMS(X2, harm);


%% Re-simulate the system for p = p+dp                              

% UPDATE SYSTEM FOR NlVib _________________________________________________

% update stiffness values
k11 = k11 + dp(1);
k13 = k13 + dp(2);
k21 = k21 + dp(3);
k23 = k23 + dp(4);
m1 = m1 + dp(5);
m2 = m2 + dp(6);
d1 = d1 + dp(7);
d2 = d2 + dp(8);

% 1st cubic spring
nonlinear_elements{1} = struct('type', 'cubicspring',...
    'stiffness', k13, 'force_direction', [1;0]);

% 2nd cubic spring
nonlinear_elements{2} = struct('type', 'cubicspring',...
    'stiffness', k23, 'force_direction', [1;-1]);

% generate system
system = ChainOfOscillators([m1,m2], [d1,d2,0], [k11,k21,0], ...
    nonlinear_elements, [F,0]);

% EIGENVALUE ANALYSIS _____________________________________________________
% stiffness matrix
K = system.K;

% mass matrix
M = system.M;

% eigenvalue analysis
[VMn,om] = eigs(K, M, 2, 'SM');
[om_n,ind] = sort(sqrt(diag(om)));

% RUN ANALYSIS_____________________________________________________________

% initial guess for continuation (linear approximation)
Q1 = (-Om_s^2*system.M + 1i*Om_s*system.D + system.K)\[F,0]';
x0 = zeros((2*H+1)*length(Q1),1);
x0(length(Q1)+(1:2*length(Q1))) = [real(Q1);-imag(Q1)]; 

% % set frequency range for analysis
% Om_s = .2*om_n(foi);	% start frequency
% Om_e = 4*om_n(foi);     % end frequency 

% re-simulate the system
[X,SolInfo,Sol] = solve_and_continue(x0,...
    @(X) HB_residual(X,system,H,N), Om_s, Om_e, ds); 

% decode re-simulated solution
Xr = nlvib_decode(X, SolInfo, Sol, 'FRF', 'HB', n, H);

% compute RMS
Xr.RMS = RMS(Xr,harm);


%% Plot results                                                     

% input of 'plot_frc1.m'. which curve to plot?
frc_data = {Xn, Xr, X1, X2};

% input of 'plot_frc1.m'. for more details look function description
items2plot.subplot = [221 222 223 224];
items2plot.object = {'mod', 'mod', 'phase', 'phase'};
items2plot.dof = [1 2 1 2];
items2plot.harm = [1 1 1 1];
% input of 'plot_frc1.m'. for more details look function description
plotOptions.legend = {'nom', 're-sim', 'lin', 'quad'};
plotOptions.color = {'k', 'b', 'g', 'm'};
plotOptions.lineWidth = {1.5 1.5 1.5 1.5};
plotOptions.lineStyle = {'-','-','-.','-.'};
plotOptions.marker = {'none', 'none', 'none', 'none'};
plotOptions.fontSize = 15;

% plot modulus of first harmonic, dof 1 and 2 and RMS of harm 0,2,3,5
figure('units','normalized','position',[.2 .2 .6 .6])
plot_FRC(frc_data, items2plot, plotOptions);

% plot coefficients A & B
items2plot.object = {'A', 'A', 'B', 'B'};
figure('units','normalized','position',[.2 .2 .6 .6])
plot_FRC(frc_data, items2plot, plotOptions);



