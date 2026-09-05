function rho=rho_a(Ta,par)
% Temperature-dependent atmospheric density

rho = par.p_a./(par.R_d*Ta);
