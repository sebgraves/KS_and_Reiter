function [v,jac] = solve_valfunc_KS(ECM,c,s,param,glob,options);


if ECM == 0
    % Solve problem using Golden Search
    B                   = menufun_KS('bounds',s,[],param,glob,options);
    obj                 = @(kp)valfunc_KS(c,s,kp,param,glob,options);     
    kp                  = goldenx(obj,B(:,1),B(:,2));
    [v1,ve,Phi_keKA]    = valfunc_KS(c,s,kp,param,glob,options);

    % Report Jacobian if necessary
    if nargout==2
        jac = [glob.Phi , - param.beta*Phi_keKA;
               -glob.Emat*glob.Phia, glob.Phi];
    end
    
elseif ECM == 1
    % Solve problem using ECM
    
    % Split coefficients
    c1  = c(1:glob.Ns);
    ce  = c(glob.Ns+1:end);

    R   = menufun_KS('R',s,[],param,glob,options);
    w   = menufun_KS('w',s,[],param,glob,options);

    V1  = funeval(c1,glob.fspace,s,[1 0 0 0]);  % Partial derivative with respect to k
    V1  = max(V1,0);                            % Ensure derivative is positive

    % Use derivative to construct k_prime and keep on the grid
    kp  = min(max(0, (1-param.tau)*w.*s(:,2) + param.b*w.*(1-s(:,2)) + R.*s(:,1) - R./V1),(1-param.tau)*w.*s(:,2) + param.b*w.*(1-s(:,2)) + R.*s(:,1));
    kp  = min(kp,glob.kmax);

    Phi_k       = splibas(glob.kgrid0,0,glob.spliorder(1),kp);
    Phi_ke      = dprod(glob.Phie,Phi_k);
    Phi_keK     = dprod(glob.PhiK,Phi_ke);
    Phi_keKA    = dprod(glob.PhiA,Phi_keK);
    
    % Compute flow payoff
    F       = menufun_KS('F',s,kp,param,glob,options);

    v1  = F + param.beta*Phi_keKA*ce;
    ve  = glob.Emat*glob.Phia*c1;

    % Report Jacobian if necessary
    if nargout==2
        jac = [glob.Phi , - param.beta*Phi_keKA;
               -glob.Emat*glob.Phia, glob.Phi];
    end
end

%% Package output
v.v1    = v1;
v.ve    = ve;
v.kp    = kp;
