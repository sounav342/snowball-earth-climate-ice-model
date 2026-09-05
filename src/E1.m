function En=E1(R,Ta,To,Ts,qa,par)
% evaporation rate from overall surface to the atmosphere
En = zeros(par.nj,1);
Ts(isnan(Ts)) = 0;

dq_To = qa - q_sat(To,par);
dq_Ts = qa - q_sat(Ts,par);
rhoa = rho_a(Ta,par);

En = -par.CD_evap .* rhoa .* ((1 - R) .* dq_To + R .* dq_Ts);
En(En < 0) = 0;
