function gzeta = g_fun(z,eta,param,glob);

gzeta = exp(param.rhoA*log(z) + eta);

% Keep on the grid
gzeta = max(min(gzeta,glob.Amax),glob.Amin);

end

