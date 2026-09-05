function Qo = Q_o(R, par)
% Qo: solar flux absorbed per unit open-ocean area.
% The EBM applies the open-ocean area fraction, (1-R), outside this helper.

[~,~,alpha_o_local,R] = surface_albedo_fractional(R,par);

SW_down = par.Tau .* par.S_col(:) .* par.Q ./ 4;
Qo = zeros(par.nj,1);
openOcean = R < 1;
Qo(openOcean) = SW_down(openOcean) .* (1 - alpha_o_local(openOcean));
end
