function [A_SS, Y_SS, C_SS, I_SS, W_SS, R_SS, L_SS, K_SS] = SS_calc(lambda,param,glob,options);

A_SS    = 1;
L_SS    = param.L;

K_SS    = lambda'*glob.s_ind(:,1);

W_SS    = (1-param.alpha)*(K_SS/L_SS)^param.alpha;
R_SS    = param.alpha*(K_SS/L_SS)^(param.alpha-1) + 1 - param.delta;

Y_SS    = K_SS^param.alpha*L_SS^(1-param.alpha);
I_SS    = param.delta*K_SS;
C_SS    = Y_SS - I_SS;
