function F = F_fun(X,Xlag,param,glob,options);

% Split the vectors
c_lag       = Xlag(1:glob.Ns);
lambda_lag  = Xlag(glob.Ns+1:glob.Ns+glob.Nsf);
logA_lag    = Xlag(glob.Ns+glob.Nsf+1);
Ec_lag      = Xlag(glob.Ns+glob.Nsf+2:end);

c       = X(1:glob.Ns);
lambda  = X(glob.Ns+1:glob.Ns+glob.Nsf);
logA    = X(glob.Ns+glob.Nsf+1);
Ec      = X(glob.Ns+glob.Nsf+2:end);

% Construct value of K
K_lag   = lambda_lag'*glob.sf(:,1);

%% Solve for RHS of first set of equations
v1   = solve_valfunc_Reiter(Ec_lag,glob.s,K_lag,logA_lag,param,glob,options);
F1   = glob.Phi*c_lag - v1.ve;

%% Solve for policy function in second set of equations
glob.Phi2   = glob.Phi2f;
glob.Phi    = glob.Phif;
glob.PnK    = glob.PnKf;
v2          = solve_valfunc_Reiter(Ec_lag,glob.sf,K_lag,logA_lag,param,glob,options);

% Create Q matrix
fspacek     = fundef({'spli',glob.kgridf,0,1});
Qk          = funbas(fspacek,v2.Kp);
Q           = dprod(glob.Qe,Qk);

%% Calculate LHS - RHS of F
F2  = lambda - Q'*lambda_lag;
F3  = logA - param.rhoA*logA_lag;
F4  = c - Ec_lag;

F   = [F1;F2;F3;F4];