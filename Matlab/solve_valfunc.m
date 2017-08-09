function [v,jac] = solve_valfunc(c,s,K,param,glob,options);

% Solve problem
B                   = menufun('bounds',s,[],K,param,glob,options);
obj                 = @(Kp)valfunc(c,s,Kp,K,param,glob,options);     
Kp                  = goldenx(obj,B(:,1),B(:,2));
[v1,ve,Phi_KpZ]     = valfunc(c,s,Kp,K,param,glob,options);


% Report Jacobian if necessary
if nargout==2
    jac = glob.Phi - param.beta*glob.PnK*Phi_KpZ;
end

% Package output
v.v1        = v1;
v.ve        = ve;
v.Kp        = Kp;
