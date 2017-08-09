function [v1,ve,Phi_KpZ] = valfunc(c,s,Kp,K,param,glob,options)

% Compute flow payoff
F       = menufun('F',s,Kp,K,param,glob,options);

% Create basis matrices for continuation value
Phi_Kp      = splibas(glob.kgrid0,0,glob.spliorder(1),Kp);
Phi_KpZ     = dprod(glob.Phi2,Phi_Kp);

% Compute value
v1      = F + param.beta*Phi_KpZ*c;
ve      = glob.PnK*(F + param.beta*Phi_KpZ*c);

