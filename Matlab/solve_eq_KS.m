function eq = solve_eq(eq_SS,param,glob,options);

%% Initialize forecast rule system
b   = [0.183640;0.869577;0.273981;-0.007685];

%% Draw T exogenous aggregate productivity values
rng(12);
sim.A_series    = zeros(options.T,1);
sim.eps_series  = normrnd(0,param.sigA,options.T,1);
sim.A_series(1) = 1;
for t=1:(options.T-1)
    sim.A_series(t+1)   = g_fun(sim.A_series(t),sim.eps_series(t),param,glob);
end

%% Krusell-Smith Algorithm
for citer = 1:options.itermaxreg
    
    options.iter_num = citer;
    
    % Create \hat K' vector
    K_p     = zeros(glob.Ns,1);
    for i = 1:glob.Ns
        K_p(i) = exp(b(1) + b(2)*log(glob.s(i,3)) + b(3)*log(glob.s(i,4)) + b(4)*log(glob.s(i,3)).*log(glob.s(i,4)));
    end
    % Keep K_p on grid
    K_p     = min(max(K_p,glob.Kmin),glob.Kmax);
    
    % Create basis matrices
    glob.PhiK   = splibas(glob.Kgrid0,0,glob.spliorder(3),K_p);
    
    % Solve Value functions, given forecast rule system
    eq      = solve_cL_KS(param,glob,options);
    glob.c  = eq.c;
    options.cresult = eq.c;
    fprintf('Solved for value functions \n');

    % Initialise cross-section using SS
    sim.mu_series           = zeros(glob.nf(1)*glob.nf(2),1,options.T+1);
    sim.mu_series(:,:,1)    = eq_SS.L;
    
    sim.K_series    = zeros(options.T+1,1);
    sim.p_series    = zeros(options.T,1);
    sim.K_series(1) = sim.mu_series(:,:,1)'*glob.s_ind(:,1);
    
    % Simulate the model
    for t = 1:options.T
        
        % Find optimal policy functions, given A and K
        st          = gridmake(glob.kgridf,glob.egrid,sim.K_series(t),sim.A_series(t));
        Ns          = size(st,1);
        
        Kprime      = exp(b(1) + b(2)*log(st(:,3)) + b(3)*log(st(:,4)) + b(4)*log(st(:,3)).*log(st(:,4)));
        glob.Phie   = splibas(glob.egrid,0,glob.spliorder(2),st(:,2));
        glob.PhiK   = splibas(glob.Kgrid0,0,glob.spliorder(3),Kprime);
        glob.PhiA   = splibas(glob.Agrid0,0,glob.spliorder(4),st(:,4));

        v           = solve_valfunc_KS(1,eq.c,st,param,glob,options);
        
        % Update discretized distribution for the next period
        kp                      = min(v.kp,max(glob.kgrid));
        Qk                      = funbas(glob.fspacek,kp);
        sim.mu_series(:,:,t+1)  = dprod(glob.Qe,Qk)'*sim.mu_series(:,:,t);
        sim.K_series(t+1)       = sim.mu_series(:,:,t+1)'*glob.s_ind(:,1);
        
        if sum(sim.mu_series(:,:,t+1)) > 1.01
            error('Error. Distribution above 1.')
        elseif sum(sim.mu_series(:,:,t+1)) < 0.99
            error('Error. Distribution less than 1.')
        end
    end
    
    %% Run regression
    
    % Drop burn period
    K_series_reg    = sim.K_series(options.Terg:end-1);
    A_series_reg    = sim.A_series(options.Terg:end);
    Kp_series_reg   = sim.K_series(options.Terg+1:end);
    
    [beta,~,~,~,stats]  = regress(log(Kp_series_reg), [ones(size(Kp_series_reg,1),1), log(K_series_reg), log(A_series_reg), log(K_series_reg).*log(A_series_reg)]);
    
    %% Update forecast coefficients
    diff    = max(abs(b - beta));
    disp(diff)
    disp(b)
    disp(beta)
    if diff < options.tolreg;
        break
    end
    
    % Dampening to construct new coeffs
    b   = options.damp*beta + (1-options.damp)*b;
    
    %% Replace some matrices
    glob.Phie = splibas(glob.egrid,0,glob.spliorder(2),glob.s(:,2));
    glob.PhiA = splibas(glob.Agrid0,0,glob.spliorder(4),glob.s(:,4));
end

eq.beta = beta;

fprintf('R-squared for forecast regression');
format long
disp(stats(1))
format short