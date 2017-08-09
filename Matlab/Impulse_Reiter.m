function IRF_Reiter = Impulse_Reiter(A,B,X_SS,param,glob,options);


%% Calculate SS values
lambda_SS   = X_SS(glob.Ns+1:glob.Ns+glob.Nsf,1,1);
[A_SS, Y_SS, C_SS, I_SS, W_SS, R_SS, L, K_SS] = SS_calc(lambda_SS,param,glob,options);

%% Draw exogenous productivity shocks
eps_series      = zeros(options.T_irf+1,1);
eps_series(1)   = param.sigA;

%% Initialise X_0 = X_SS
X_diff      = zeros(size(A,1),1,options.T_irf+1);
Xt          = zeros(size(A,1),1,options.T_irf+1);
Xt(:,:,1)   = X_SS;

%% Simulate X_t
for i = 1:options.T_irf
    X_diff(:,:,i+1) = A*X_diff(:,:,i) + B*eps_series(i);
    Xt(:,:,i+1)     = X_SS + X_diff(:,:,i+1);
end

%% Calculate Aggregates
A_series    = Xt(glob.Ns+glob.Nsf+1,1,:);
A_series    = reshape(A_series,options.T_irf+1,1);
A_series    = exp(A_series);
K_series    = zeros(options.T_irf+1,1);
C_series    = zeros(options.T_irf,1);
Y_series    = zeros(options.T_irf,1);
I_series    = zeros(options.T_irf,1);
W_series    = zeros(options.T_irf,1);
R_series    = zeros(options.T_irf,1);

% Initialise K
lambda_0      = Xt(glob.Ns+1:glob.Ns+glob.Nsf,1,1);
K_series(1)   = lambda_0'*glob.sf(:,1);

for t = 1:options.T_irf
    lambda_t        = Xt(glob.Ns+1:glob.Ns+glob.Nsf,1,t+1);
    K_series(t+1)   = lambda_t'*glob.sf(:,1);
    
    Y_series(t)     = A_series(t)*K_series(t)^param.alpha*L^(1-param.alpha);
    I_series(t)     = K_series(t+1) - (1-param.delta)*K_series(t);
    C_series(t)     = Y_series(t) - I_series(t);
    W_series(t)     = (1-param.alpha)*A_series(t).*(K_series(t)/L).^param.alpha;
    R_series(t)     = param.alpha*A_series(t).*(K_series(t)/L).^(param.alpha-1) + 1 - param.delta;
end


%% Save IRFs

IRF_Reiter.A   = A_series(2:end) - A_SS;
IRF_Reiter.C   = log(C_series(2:end)) - log(C_SS);
IRF_Reiter.Y   = log(Y_series(2:end)) - log(Y_SS);
IRF_Reiter.I   = log(I_series(2:end)) - log(I_SS);
IRF_Reiter.W   = log(W_series(2:end)) - log(W_SS);
IRF_Reiter.R   = R_series(2:end) - R_SS;
