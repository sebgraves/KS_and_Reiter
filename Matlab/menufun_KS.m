function out = menufun(flag,s,Kp,param,glob,options)

% Parameters
alpha   = param.alpha;
beta    = param.beta;
delta   = param.delta;
b       = param.b;
L       = param.L;
tau     = param.tau;

k   = s(:,1);
e   = s(:,2);

% Equilibrium prices
w           = (1-alpha)*s(:,4).*(s(:,3)/L).^alpha;
R           = alpha*s(:,4).*(s(:,3)/L).^(alpha-1) + 1 - delta;

switch flag
    case 'bounds'
        kpmin   = ones(size(s,1),1)*min(glob.kgrid);
        kpmax   = w*(1-tau).*e + b*w.*(1-e) + R.*k;
        kpmax   = min(glob.kmax,kpmax);
        out     = [kpmin,kpmax];
    case 'F'
        out = log(w*(1-tau).*e + b*w.*(1-e) + R.*k - Kp);
    case 'Cons'
        out = w*(1-tau).*e + b*w.*(1-e) + R.*k - Kp;
    case 'w'
        out = w;
    case 'R'
        out = R;
end

end
