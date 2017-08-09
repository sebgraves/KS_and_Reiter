function [v,jac] = solve_valfunc(c,s,K,logA,param,glob,options);

% Solve problem
B                   = menufun_Reiter('bounds',s,[],K,logA,param,glob,options);
obj                 = @(Kp)valfunc_Reiter(c,s,Kp,K,logA,param,glob,options);     
Kp                  = goldenx(obj,B(:,1),B(:,2));
[v1,ve,Phi_KpZ]     = valfunc_Reiter(c,s,Kp,K,logA,param,glob,options);


% Report Jacobian if necessary
if nargout==2
    jac = glob.Phi - param.beta*glob.PnK*Phi_KpZ;
end

% Package output
v.v1        = v1;
v.ve        = ve;
v.Kp        = Kp;
