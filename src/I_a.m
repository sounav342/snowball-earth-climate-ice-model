function Ia=I_a(R,par,Ta,To,Ts,df,n)
% Ia: grid-cell mean thermal infrared flux absorbed by the atmosphere.

R = validate_fractional_ice_cover(R,par,'I_a');
Ta = Ta(:);
To = To(:);
Ts = Ts(:);
Ts(isnan(Ts)) = 0;

eps = par.A_a;
% Area-weight upward longwave emission from open ocean and ice fractions.
Ia = (1 - R) .* eps .* par.sigma .* To.^4 ...
    + R .* eps .* par.sigma .* Ts.^4 ...
    - 2 .* eps .* par.sigma .* Ta.^4;
end
