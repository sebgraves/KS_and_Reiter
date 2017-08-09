function Sim_Reiter = Simulate(A,B,X_SS,param,glob,options);

%% Draw exogenous productivity shocks
rng(100)
eps_series  = normrnd(0,param.sigA,options.T,1);

%% Initialise X_0 = X_SS
X_diff      = zeros(size(A,1),1,options.T+1);
Xt          = zeros(size(A,1),1,options.T+1);
Xt(:,:,1)   = X_SS;

%% Simulate X_t
for i = 1:options.T
    X_diff(:,:,i+1) = A*X_diff(:,:,i) + B*eps_series(i);
    Xt(:,:,i+1)     = X_SS + X_diff(:,:,i+1);
end

%% Calculate Aggregates
A_series    = Xt(glob.Ns+glob.Nsf+1,1,:);
A_series    = reshape(A_series,options.T+1,1);
A_series    = exp(A_series);
K_series    = zeros(options.T+1,1);
C_series    = zeros(options.T,1);
Y_series    = zeros(options.T,1);
I_series    = zeros(options.T,1);
W_series    = zeros(options.T,1);
R_series    = zeros(options.T,1);

% Initialise K
lambda_0      = Xt(glob.Ns+1:glob.Ns+glob.Nsf,1,1);
K_series(1)   = lambda_0'*glob.sf(:,1);

for t = 1:options.T
    lambda_t        = Xt(glob.Ns+1:glob.Ns+glob.Nsf,1,t+1);
    K_series(t+1)   = lambda_t'*glob.sf(:,1);
    
    Y_series(t)     = A_series(t)*K_series(t)^param.alpha*param.L^(1-param.alpha);
    I_series(t)     = K_series(t+1) - (1-param.delta)*K_series(t);
    C_series(t)     = Y_series(t) - I_series(t);
    W_series(t)     = (1-param.alpha)*A_series(t).*(K_series(t)/param.L).^param.alpha;
    R_series(t)     = param.alpha*A_series(t).*(K_series(t)/param.L).^(param.alpha-1) + 1 - param.delta;
end

%% Save Simulations (dropping burn period)
Sim_Reiter.A   = A_series(options.Terg:end);
Sim_Reiter.C   = C_series(options.Terg:end);
Sim_Reiter.Y   = Y_series(options.Terg:end);
Sim_Reiter.I   = I_series(options.Terg:end);
Sim_Reiter.W   = W_series(options.Terg:end);
Sim_Reiter.R   = R_series(options.Terg:end);
Sim_Reiter.K   = K_series(options.Terg:end);

%% Make Table Output

[~,A] = hpfilter(log(Sim_Reiter.A),100);
[~,C] = hpfilter(log(Sim_Reiter.C),100);
[~,Y] = hpfilter(log(Sim_Reiter.Y),100);
[~,I] = hpfilter(log(Sim_Reiter.I),100);
[~,W] = hpfilter(log(Sim_Reiter.W),100);
[~,R] = hpfilter(Sim_Reiter.R,100);

% Make Table 2
Col1        = [std(Y);std(C)/std(Y);std(I)/std(Y);std(W)/std(Y);std(R)/std(Y)];
CorrMatrix  = corrcoef([Y,C,I,W,R]);
Col2        = [CorrMatrix(1,1);CorrMatrix(1,2);CorrMatrix(1,3);CorrMatrix(1,4);CorrMatrix(1,5)];

Sim_Reiter.TableData   = round([Col1,Col2],4);
