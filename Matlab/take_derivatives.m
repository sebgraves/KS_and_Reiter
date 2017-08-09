function [F1,F2,F3,F4] = take_derivatives(X_SS,param,glob,options);

Nx      = size(X_SS,1);

%% Exact calculation of F3 and F4
F3                  = zeros(Nx,1);
F3(end - glob.Ns)   = -1;

F4                  = [zeros(Nx - glob.Ns,glob.Ns);-eye(glob.Ns)];

%% Approximate F1 and F2
F1      = zeros(Nx,Nx);
F2      = zeros(Nx,Nx);
F_SS    = F_fun(X_SS,X_SS,param,glob,options);
deriv   = options.reltol;

for i = 1:Nx
    Xu  = X_SS;
	h = X_SS(i) * deriv;
	if (h < deriv)
		h = deriv;
	end
	Xu(i) = Xu(i) + h;
	Fu1 = F_fun(Xu,X_SS,param,glob,options);
    Fu2 = F_fun(X_SS,Xu,param,glob,options);
	F1(:,i) = (Fu1 - F_SS)./h;
    F2(:,i) = (Fu2 - F_SS)./h;
end


    
    