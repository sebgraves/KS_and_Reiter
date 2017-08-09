function [v1,ve,Phi_keKA] = valfunc(c,s,kp,param,glob,options)

% Split coefficients for c1,ce
c1  = c(1:glob.Ns);
ce  = c(glob.Ns+1:end);

% Compute flow payoff
F       = menufun_KS('F',s,kp,param,glob,options);

% Create basis matrices for continuation value
Phi_k       = splibas(glob.kgrid0,0,glob.spliorder(1),kp);
Phi_ke      = dprod(glob.Phie,Phi_k);
Phi_keK     = dprod(glob.PhiK,Phi_ke);
Phi_keKA    = dprod(glob.PhiA,Phi_keK);

% Compute values
v1      = F + param.beta*Phi_keKA*ce;
ve      = glob.Emat*glob.Phia*c1;
