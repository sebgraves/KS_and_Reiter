%% Solve Incomplete Markets Model using Krusell-Smith and Reiter Methods

%% Set Options
% All
options.Nbell       = 5;        % Number of Bellman iterations before Newton.
options.Nnewt       = 30;       % Maximum number of Newton steps
options.tolc        = 1e-8;     % Tol on value functions
options.T_irf       = 50;       % Number of periods for IRFs
options.T           = 1000;     % Number of periods for simulation
options.Terg        = 500;      % Number of periods to burn

% Stationary
options.itermaxL    = 5000;     % Max number of iterations to find L
options.tolL        = 1e-11;    % Tol for L
options.tolK        = 1e-5;     % Tol for equilibrium K
options.itermaxK    = 100;      % Max iterations for bisection
% Reiter
options.reltol      = 1e-6;     % Relative tol for derivatives
% KS
options.tolreg      = 1e-6;     % Tol for forecast coefficients
options.itermaxreg  = 200;      % Max number of iterations to find coeffs
options.damp        = 0.4;      % Weight on **new** forecast rule coefficients in fixed-point iteration
options.cresult     = [];

%% Set globals and parameters
glob.n          = [40,2,6,6];   % Number of nodes for k e, K, A
glob.nf         = [100,2];      % Number of nodes for distribution histogram
glob.spliorder  = [2,1,1,1];    % Order of splines (Envelope Condition Method seems only robust when quadratic in k)
glob.Na         = 30;           % Number of nodes for quadrature
glob.curv       = 1;            % Amount of curvature for k grid
glob.kmin       = 0;            % Lower bound for k
glob.kmax       = 15;           % Upper bound for k
glob.Kmin       = 3.4;          % Lower bound for K
glob.Kmax       = 4.6;          % Upper bound for K
glob.Amin       = 0.9;          % Lower bound for A (chosen based on sim, should use better rule based on stationary dist?)
glob.Amax       = 1.1;          % Upper bound for A (chosen based on sim)

%% Model parameters
param.alpha     = 0.36;
param.beta      = 0.96;
param.delta     = 0.1;
param.b         = 0.15;
param.rhoA      = 0.859;
param.sigA      = 0.014;
param.pi01      = 0.5;
param.pi10      = 0.038;
param.L         = 1 - (param.pi10)/(param.pi01 + param.pi10);
param.tau       = param.b*(1-param.L)/param.L;

%% Solve for Stationary Equilibrium
[param,glob]    = setup(param,glob,options);
eq_SS           = solve_eq(param,glob,options);
fprintf('Solved stationary equilibrium\n');

%% Solve Reiter
eq_Reiter       = solve_eq_Reiter(eq_SS,param,glob,options);
fprintf('Solved full model using Reiter method \n');
IRF_Reiter      = Impulse_Reiter(eq_Reiter.A,eq_Reiter.B,eq_Reiter.X_SS,param,glob,options);
Sim_Reiter      = Simulate_Reiter(eq_Reiter.A,eq_Reiter.B,eq_Reiter.X_SS,param,glob,options);
fprintf('Calculated IRFs and Simulated using Reiter method \n');

%% Solve KS
[param,glob]    = setup_KS(param,glob,options);   
eq_KS           = solve_eq_KS(eq_SS,param,glob,options);
fprintf('Solved full model using KS method \n');
IRF_KS          = Impulse_KS(eq_KS,eq_SS,param,glob,options);
Sim_KS          = Simulate_KS(eq_KS,eq_SS,param,glob,options);
fprintf('Calculated IRFs and Simulated using KS method \n');

%% Replicate Figure 1 and 2 From Winberry (2016) User Guide

% Replicate Figures 1 and 2 from Winberry User Guide
M = figure(105);
set(M,'Pos',[100          100        400         600]);
subplot(311);
plot(glob.kgridf,eq_SS.Cons_e,glob.kgridf,eq_SS.Cons_u)
title('Consumption Decision Rule')
xlabel('Assets, a')
ylabel('Consumption, c(a,e)')
legend('Employed','Unemployed')
legend('Location','southeast')
subplot(312);
plot(glob.kgridf,eq_SS.Kp_e,glob.kgridf,eq_SS.Kp_u)
title('Saving Decision Rule')
xlabel('Assets, a')
ylabel('Saving, g(a,e)')
legend('Employed','Unemployed')
legend('Location','southeast')
subplot(313);
hold on;
plot(glob.kgridf,eq_SS.L_e)
ylabel('Mass of households')
yyaxis right
plot(glob.kgridf,eq_SS.L_u)
yyaxis right
title('Invariant Distribution of Households')
xlabel('Assets, a')
legend('Employed (LHS)','Unemployed (RHS)')
legend('Location','southeast')
print('SS_Plots','-djpeg')

%% Replicate Figure 3 from Winberry User Guide

H = figure(111);
set(H,'Pos',[100 100 600 800]);
subplot(321);
hold on; plot(IRF_KS.A,'r','LineWidth',1); plot(IRF_Reiter.A,'b','LineWidth',1)
legend({'KS','Reiter'})
title('Aggregate TFP')
subplot(322);
hold on; plot(IRF_KS.Y,'r','LineWidth',1); plot(IRF_Reiter.Y,'b','LineWidth',1)
title('log Y')
subplot(323);
hold on; plot(IRF_KS.C,'r','LineWidth',1); plot(IRF_Reiter.C,'b','LineWidth',1)
title('log C')
subplot(324);
hold on; plot(IRF_KS.I,'r','LineWidth',1); plot(IRF_Reiter.I,'b','LineWidth',1)
title('log I')
subplot(325);
hold on; plot(IRF_KS.W,'r','LineWidth',1); plot(IRF_Reiter.W,'b','LineWidth',1)
title('log W')
subplot(326);
hold on; plot(IRF_KS.R,'r','LineWidth',1); plot(IRF_Reiter.R,'b','LineWidth',1)
title('R')
print('IRF_Compare','-djpeg')

% Plot difference in IRFs as they are so close
G = figure(222);
set(G,'Pos',[100 100 400 300]);
hold on; 
plot(IRF_KS.Y-IRF_Reiter.Y,'r','LineWidth',1);
plot(IRF_KS.C-IRF_Reiter.C,'b','LineWidth',1);
plot(IRF_KS.I-IRF_Reiter.I,'k','LineWidth',1);
plot(IRF_KS.R-IRF_Reiter.R,'g','LineWidth',1);
legend({'Y','C','I','W','R'})
print('IRF_Diffs','-djpeg')

%% Plot representative sims
J = figure(234);
set(J,'Pos',[100 100 400 300]);
hold on; 
plot(Sim_KS.K,'r','LineWidth',1);
plot(Sim_Reiter.K,'b','LineWidth',1);
legend({'KS','Reiter'})
xlabel('t')
ylabel("K")
print('Sim_KS_Reiter','-djpeg')



