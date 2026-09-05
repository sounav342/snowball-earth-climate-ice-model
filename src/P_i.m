function Pi=P_i(Pn,Ts,par)
% Rate of Snowfall in prognostic ice h eq in kg/m^2/s
Pi = zeros(par.nj,1);
index = ~isnan(Ts);
Pi(index) = Pn(index);
