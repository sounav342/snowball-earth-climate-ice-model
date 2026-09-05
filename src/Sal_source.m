function [Sal_Sn,dhdt_cond,dhdt_odiff] = Sal_source(dhdt_odiff,h,R,Ts,T_f,par)
% This subroutine ensures snowfall occurs over ice and rain over ocean,
% in the source function submitted to the ice flow model.
% Conductive Flux

C = Cflux(h,R,Ts,T_f,par);
dhdt_cond = C/(par.rho_o*par.L_f);

Sal_Sn = zeros(par.nj,1);
Sal_Sn = dhdt_cond + dhdt_odiff;

% Boundary Conditions on Source Term

end
