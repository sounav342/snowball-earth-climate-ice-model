function Hs=H_s(R,Ta,Ts,par)
% Sensible heat flux per unit ice-covered area.
% The ice surface equation uses this per-ice-area flux where R > 0.

R = validate_fractional_ice_cover(R,par,'H_s');
Ta = Ta(:);
Ts = Ts(:);
Ts(isnan(Ts)) = 0;
rhoa = rho_a(Ta,par);
Hs = zeros(par.nj,1);

iceCovered = R > 0;
Hs(iceCovered) = -rhoa(iceCovered) .* par.c_a .* par.CD_sens ...
    .* (Ts(iceCovered) - Ta(iceCovered));
end
