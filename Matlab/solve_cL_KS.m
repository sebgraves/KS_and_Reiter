function eq = solve_cL(param,glob,options);

%% Initialise guess using constant assets

if isempty(options.cresult)==0;
    cold = options.cresult;
else
    c1old   = zeros(glob.Ns,1);
    ceold   = zeros(glob.Ns,1);
    cold    = [c1old;ceold];
end

%% For first round, solve using Golden Search
if options.iter_num == 1

    %% Bellman iterations, Golden Search
    for citer = 1:options.Nbell
        % 1. Compute values
        v   = solve_valfunc_KS(0,cold,glob.s,param,glob,options);
        % 2. Update c
        c1  = glob.Phi\v.v1;
        ce  = glob.Phi\v.ve;
        c   = [c1;ce];
        % 3. Compute distance and update
        dc      = norm(c-cold)/norm(cold);
        cold    = c;
        if dc < options.tolc
            break
        end
    end

    %% Newton iterations, Golden Search
    for citer = 1:(options.Nnewt)
        % 1. Compute values
        [v,jac]     = solve_valfunc_KS(0,cold,glob.s,param,glob,options);
        % 2. Update c
        c1old       = cold(1:glob.Ns);
        ceold       = cold(glob.Ns+1:end);
        c           = cold - jac\[glob.Phi*c1old - v.v1;
                                  glob.Phi*ceold - v.ve];
        % 3. Compute distances and update
        dc      = norm(c-cold)/norm(cold);
        cold    = c;
        % 4. Check convergence
        if dc < options.tolc
            break
        end
        if citer == options.Nnewt
            fprintf('Newton iteration failed \n');
            return
        end
    end
else
    %% For subsequent rounds, use Newton iterations with Envelope Condition Method
    for citer = 1:(options.Nnewt)
        % 1. Compute values
        [v,jac]     = solve_valfunc_KS(1,cold,glob.s,param,glob,options);
        % 2. Update c
        c1old       = cold(1:glob.Ns);
        ceold       = cold(glob.Ns+1:end);
        c           = cold - jac\[glob.Phi*c1old - v.v1;
                                  glob.Phi*ceold - v.ve];
        % 3. Compute distances and update
        dc      = norm(c-cold)/norm(cold);
        cold    = c;
        % 4. Check convergence
        if dc < options.tolc
            break
        end
        if citer == options.Nnewt
            fprintf('Newton iteration failed \n');
            return
        end
    end
end

%% Save solutions
options.cresult = c;

%% Pack up output
eq.c    = c;
eq.v    = v;
end