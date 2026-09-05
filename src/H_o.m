function Ho=H_o(R,Ta,To,par)
% Sensible heat flux per unit open-ocean area.
% The EBM applies the open-ocean area fraction, (1-R), outside this helper.

R = validate_fractional_ice_cover(R,par,'H_o');
Ta = Ta(:);
To = To(:);
rhoa = rho_a(Ta,par);
Ho = zeros(par.nj,1);

openOcean = R < 1;
Ho(openOcean) = -rhoa(openOcean) .* par.c_a .* par.CD_sens ...
    .* (To(openOcean) - Ta(openOcean));
end
