function [v1,ve,Phi_KpZ] = valfunc(c,s,Kp,K,logA,param,glob,options)

% Compute flow payoff
F       = menufun_Reiter('F',s,Kp,K,logA,param,glob,options);

% Create basis matrices for continuation value
Phi_Kp      = splibas(glob.kgrid0,0,glob.spliorder(1),Kp);
Phi_KpZ     = dprod(glob.Phi2,Phi_Kp);

v1      = F + param.beta*Phi_KpZ*c;
ve      = glob.PnK*(F + param.beta*Phi_KpZ*c);

