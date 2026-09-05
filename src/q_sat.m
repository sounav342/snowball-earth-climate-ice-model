function qsat=q_sat(Ta,par)
% saturation specific humidity over water
% Pa in mb, T in K
% % first saturation mixing ratio over water (gr water vapor per gram dry air):

rw = r_w(Ta,par);
% and from that, saturation specific humidity (gr water vapor per gram moist air):
qsat = rw./(1+rw);
