function Sim_KS = Simulate_KS(eq_KS,eq_SS,param,glob,options);

%% Import forecast rule coefficients
b   = eq_KS.beta;

%% Draw exogenous aggregate productivity values
rng(100)
sim.A_series    = zeros(options.T,1);
sim.eps_series  = normrnd(0,param.sigA,options.T,1);
sim.A_series(1) = 1;
for t=1:(options.T-1)
    sim.A_series(t+1)   = exp(param.rhoA*log(sim.A_series(t)) + sim.eps_series(t));
end

%% Simulate endogenous variables

% Initialise cross-section using uniform distribution
sim.mu_series           = zeros(glob.nf(1)*glob.nf(2),1,options.T+1);
sim.mu_series(:,:,1)    = eq_SS.L;
    
sim.K_series    = zeros(options.T+1,1);
sim.K_series(1) = sim.mu_series(:,:,1)'*glob.s_ind(:,1);
sim.C_series    = zeros(options.T,1);
sim.Y_series    = zeros(options.T,1);
sim.I_series    = zeros(options.T,1);
sim.W_series    = zeros(options.T,1);
sim.R_series    = zeros(options.T,1);
    
% Simulate the model
for t = 1:options.T
    % Find optimal policy functions, given A and K
    st          = gridmake(glob.kgridf,glob.egrid,sim.K_series(t),sim.A_series(t));
    Ns          = size(st,1);
    Kprime      = exp(b(1) + b(2)*log(st(:,3)) + b(3)*log(st(:,4)) + b(4)*log(st(:,3)).*log(st(:,4)));
    glob.Phie   = splibas(glob.egrid,0,glob.spliorder(2),st(:,2));
    glob.PhiK   = splibas(glob.Kgrid0,0,glob.spliorder(3),Kprime);
    glob.PhiA   = splibas(glob.Agrid0,0,glob.spliorder(4),st(:,4));

    v           = solve_valfunc_KS(1,eq_KS.c,st,param,glob,options);
        
    % Update discretized distribution for the next period
    kp                      = min(v.kp,max(glob.kgrid));
    Qk                      = funbas(glob.fspacek,kp);
    sim.mu_series(:,:,t+1)  = dprod(glob.Qe,Qk)'*sim.mu_series(:,:,t);
    sim.K_series(t+1)       = sim.mu_series(:,:,t+1)'*glob.s_ind(:,1);
    
    % Update other variables
    sim.Y_series(t)     = sim.A_series(t)*sim.K_series(t)^param.alpha*param.L^(1-param.alpha);
    sim.I_series(t)     = sim.K_series(t+1) - (1-param.delta)*sim.K_series(t);
    sim.C_series(t)     = sim.Y_series(t) - sim.I_series(t);
    sim.W_series(t)     = (1-param.alpha)*sim.A_series(t).*(sim.K_series(t)/param.L).^param.alpha;
    sim.R_series(t)     = param.alpha*sim.A_series(t).*(sim.K_series(t)/param.L).^(param.alpha-1) + 1 - param.delta;
        
    if sum(sim.mu_series(:,:,t+1)) > 1.01
        error('Error. Distribution above 1.')
    elseif sum(sim.mu_series(:,:,t+1)) < 0.99
        error('Error. Distribution less than 1.')
    end
end

%% Save Simulations (dropping burn period)
Sim_KS.A   = sim.A_series(options.Terg:end);
Sim_KS.C   = sim.C_series(options.Terg:end);
Sim_KS.Y   = sim.Y_series(options.Terg:end);
Sim_KS.I   = sim.I_series(options.Terg:end);
Sim_KS.W   = sim.W_series(options.Terg:end);
Sim_KS.R   = sim.R_series(options.Terg:end);
Sim_KS.K   = sim.K_series(options.Terg:end);

%% Make Table Output

[~,A] = hpfilter(log(Sim_KS.A),100);
[~,C] = hpfilter(log(Sim_KS.C),100);
[~,Y] = hpfilter(log(Sim_KS.Y),100);
[~,I] = hpfilter(log(Sim_KS.I),100);
[~,W] = hpfilter(log(Sim_KS.W),100);
[~,R] = hpfilter(Sim_KS.R,100);

% Make Table 2
Col1        = [std(Y);std(C)/std(Y);std(I)/std(Y);std(W)/std(Y);std(R)/std(Y)];
CorrMatrix  = corrcoef([Y,C,I,W,R]);
Col2        = [CorrMatrix(1,1);CorrMatrix(1,2);CorrMatrix(1,3);CorrMatrix(1,4);CorrMatrix(1,5)];

Sim_KS.TableData   = round([Col1,Col2],4);