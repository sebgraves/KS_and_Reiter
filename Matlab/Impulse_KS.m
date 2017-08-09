function IRF_KS = Impulse_KS(eq_KS,eq_SS,param,glob,options);

%% Import forecast rule coefficients
b   = eq_KS.beta;

%% Calculate SS values
[A_SS, Y_SS, C_SS, I_SS, W_SS, R_SS, L] = SS_calc(eq_SS.L,param,glob,options);

%% Draw exogenous aggregate productivity values
sim.A_series    = ones(options.T_irf+101,1);
sim.A_series(102) = exp(param.sigA);
for t=103:options.T_irf+101
    sim.A_series(t) = exp(param.rhoA*log(sim.A_series(t-1)));
end

%% Simulate endogenous variables

% Initialise cross-section using SS
sim.mu_series           = zeros(glob.nf(1)*glob.nf(2),1,options.T+1);
sim.mu_series(:,:,1)    = eq_SS.L;
    
sim.K_series    = zeros(options.T_irf+101,1);
sim.K_series(1) = sim.mu_series(:,:,1)'*glob.s_ind(:,1);
sim.C_series    = zeros(options.T_irf+100,1);
sim.Y_series    = zeros(options.T_irf+100,1);
sim.I_series    = zeros(options.T_irf+100,1);
sim.W_series    = zeros(options.T_irf+100,1);
sim.R_series    = zeros(options.T_irf+100,1);
    
% Simulate the model
for t = 1:options.T_irf+100
    % Find optimal policy functions, given A and K
    st          = gridmake(glob.kgridf,glob.egrid,sim.K_series(t),sim.A_series(t));
    Ns          = size(st,1);
    Kprime      = exp(b(1) + b(2)*log(sim.K_series(t)) + b(3)*log(sim.A_series(t)) + b(4)*log(sim.K_series(t)).*log(sim.A_series(t)));
    glob.Phie   = splibas(glob.egrid,0,glob.spliorder(2),st(:,2));
    glob.PhiK   = splibas(glob.Kgrid0,0,glob.spliorder(3),Kprime*ones(size(st,1),1));
  	glob.PhiA   = splibas(glob.Agrid0,0,glob.spliorder(4),sim.A_series(t)*ones(size(st,1),1));

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


%% Save IRFs

logC = log(sim.C_series);
logY = log(sim.Y_series);
logI = log(sim.I_series);
logW = log(sim.W_series);

IRF_KS.A = sim.A_series(102:end) - sim.A_series(101);
IRF_KS.C = logC(102:end) - logC(101);
IRF_KS.Y = logY(102:end) - logY(101);
IRF_KS.I = logI(102:end) - logI(101);
IRF_KS.W = logW(102:end) - logW(101);
IRF_KS.R = sim.R_series(102:end) - sim.R_series(101);
        