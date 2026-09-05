function dhdt=Cforc(h,R,Ts,T_f,par)
% Conductive Basal Forcing of Sea Ice Melt/Freezing (used in ice flow model)
% The output90-par.theta(je:-1:jb),S(je:-1:jb)*par.year*100,'--ob', of this function has units of m/s.

C = Cflux(h,R,Ts,T_f,par);
dhdt = C/(par.rho_i*par.L_f);
