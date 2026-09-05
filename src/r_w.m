function rw=r_w(Ta,par)
% saturation mixing ratio over water
% Pa in mb, bc so is ew(T_a)
% T in kelvin, as that's what required by ew(T_a)
% output mixing ratio is gr/gr.
ew = e_w(Ta,par);

par.epsilon = par.R_d/par.R_v;

P = Pa_2_mb(par.p_a);
rw = par.epsilon*ew./(P-ew);
