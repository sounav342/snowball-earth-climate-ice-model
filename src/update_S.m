function [Sn,dhdt_cond] = update_S(dhdt_odiff,h,R,S,Ts,T_f,par)

% This subroutine ensures snowfall occurs over ice and rain over ocean,
% in the source function submitted to the ice flow model.
% Conductive Flux

dhdt_cond = Cforc(h,R,Ts,T_f,par);

Sn = zeros(par.nj,1);
Sn(~isnan(Ts)) = S(~isnan(Ts));

Sn = dhdt_cond + dhdt_odiff;

% Boundary Conditions on Source Term
Sn(1) = 0;
Sn(end) = 0;

end
