function Ha=H_a(R,Ho,Hs,par)
% Grid-cell mean sensible heat flux from the overall surface to atmosphere.

R = validate_fractional_ice_cover(R,par,'H_a');
Ho = Ho(:);
Hs = Hs(:);
% Area-weight ocean and ice sensible heat fluxes over the fractional cell.
Ha = -((1 - R) .* Ho + R .* Hs);
end
