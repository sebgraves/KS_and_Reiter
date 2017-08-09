function eq = solve_eq(param,glob,options)

Klb    = glob.Kmin;
Kub    = glob.Kmax;

% Storage
Kinvec  = zeros(options.itermaxK,1);
Koutvec = zeros(options.itermaxK,1);
options.cresult = [];
tictic = tic;
for tt = 1:options.itermaxK
    % 1. Update p
    K   = 0.5*(Klb+Kub);
    % 2. Solve economy, given K
    eq  = solve_cL(K,param,glob,options);
    options.cresult = eq.c;     % Save to use as starting guess in future
    % 3. Record output and print
    Kinvec(tt)  = K;
    Koutvec(tt) = eq.K;
    %fprintf('%2i. Kin:\t%2.6f\t Kout:\t%2.6f\tt:%2.1f\n',tt,K,eq.K,toc(tictic));
    % 4. Set some flags
    eq.flag.equi    = abs(Kinvec(tt)-Koutvec(tt))<options.tolK;
    eq.flag.down    = Kinvec(tt)>Koutvec(tt);
    eq.flag.up      = Kinvec(tt)<Koutvec(tt);
    % 5. Shift bounds
    Klb     = eq.flag.up*K + eq.flag.down*Klb;
    Kub     = eq.flag.up*Kub + eq.flag.down*K;
    % Break if equilibrium
    if eq.flag.equi
        break
    end
end

%% Replicate fig 1 and 2 from Winberry (2016) User Guide
eq.Kp_u = eq.v.Kp(1:glob.Nkf);
eq.Kp_e = eq.v.Kp(glob.Nkf+1:end);

Cons   = menufun('Cons',glob.sf,eq.v.Kp,eq.K,param,glob,options);
eq.Cons_u = Cons(1:glob.Nkf);
eq.Cons_e = Cons(glob.Nkf+1:end);

eq.L_u = eq.L(1:glob.Nkf);
eq.L_e = eq.L(glob.Nkf+1:end);