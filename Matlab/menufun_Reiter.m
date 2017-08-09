function out = menufun(flag,s,Kp,K,logA,param,glob,options)

% Parameters
alpha   = param.alpha;
beta    = param.beta;
delta   = param.delta;
b       = param.b;

% Equilibrium variables
L           = 1 - (param.pi10)/(param.pi01 + param.pi10);
tau         = param.b*(1-L)/L;
w           = (1-alpha)*exp(logA)*(K/L)^alpha;
R           = alpha*exp(logA)*(K/L)^(alpha-1) + 1 - delta;

switch flag
    case 'bounds'
        k   = s(:,1);
        e   = s(:,2);
        kpmin   = ones(size(s,1),1)*min(glob.kgrid);
        kpmax   = w*(1-tau)*e + b*w*(1-e) + R*k;
        kpmax   = min(glob.kmax,kpmax);
        out     = [kpmin,kpmax];
    case 'F'
        k   = s(:,1);
        e   = s(:,2);
        out = log(w*(1-tau)*e + b*w*(1-e) + R*k - Kp);
    case 'Cons'
        k   = s(:,1);
        e   = s(:,2);
        out = w*(1-tau)*e + b*w*(1-e) + R*k - Kp;
end

end
