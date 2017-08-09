function eq = solve_cL(K,param,glob,options);

%% Initialise guess
cold   = zeros(glob.Ns,1);

if isempty(options.cresult)==0;
    cold = options.cresult;
end

totaltic    = tic;

%% Bellman iterations
%fprintf('~~~~~ Bellman iterations ~~~~~\n');
for citer = 1:options.Nbell
    % 1. Compute values
    v   = solve_valfunc(cold,glob.s,K,param,glob,options);
    % 2. Update c
    c  = glob.Phi\v.ve;
    % 3. Compute distance and update
    dc      = norm(c-cold)/norm(cold);
    cold    = c;
end

%% Newton iterations
%fprintf('~~~~~ Newton iterations ~~~~~\n');
for citer = 1:(options.Nnewt)
    % 1. Compute values
    [v,jac]     = solve_valfunc(cold,glob.s,K,param,glob,options);
    % 2. Update c
    c           = cold - jac\[glob.Phi*cold - v.ve];
    % 3. Compute distances and update
    dc      = norm(c-cold)/norm(cold);
    cold    = c;
    % 4. Check convergence
    if dc < options.tolc
        break
    end
end

%% Solve again on a finer grid for K
glob.Phi2   = glob.Phi2f;
glob.Phi    = glob.Phif;
glob.PnK    = glob.PnKf;
v           = solve_valfunc(c,glob.sf,K,param,glob,options);

%% Compute stationary distribution
fspacek     = fundef({'spli',glob.kgridf,0,1});
Qk          = funbas(fspacek,v.Kp);
Q           = dprod(glob.Qe,Qk);

% Initialise uniform histogram
L           = ones(size(Q,1),1);
L           = L/sum(L);

% Iterate to convergence
for itL = 1:options.itermaxL;
    Lnew    = Q'*L;
    dL      = norm(Lnew - L)/norm(L);
    if dL<options.tolL
        break
    end
    L   = Lnew;
end

%% Compute implied K
K1      = L'*v.Kp;

%% Pack up output
eq.v    = v;
eq.c    = c;
eq.L    = L;
eq.Q    = Q;
eq.Qk   = Qk;
eq.K    = K1;
end

