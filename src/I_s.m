function Is=I_s(R,Ta,Ts,df,par)
% Is: net thermal infrared flux absorbed per unit ice-covered area.
% The ice surface equation uses this per-ice-area flux where R > 0.

R = validate_fractional_ice_cover(R,par,'I_s');
Ta = Ta(:);
Ts = Ts(:);
Ts(isnan(Ts))=0;
eps = par.A_a;
Is = zeros(par.nj,1);

iceCovered = R > 0;
Is(iceCovered) = eps .* par.sigma .* Ta(iceCovered).^4 ...
    - par.sigma .* Ts(iceCovered).^4;
end
