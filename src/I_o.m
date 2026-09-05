function Io=I_o(R,Ta,To,df,par)
% Io: net thermal infrared flux absorbed per unit open-ocean area.
% The EBM applies the open-ocean area fraction, (1-R), outside this helper.

R = validate_fractional_ice_cover(R,par,'I_o');
Ta = Ta(:);
To = To(:);
eps = par.A_a;
Io = zeros(par.nj,1);

openOcean = R < 1;
Io(openOcean) = eps .* par.sigma .* Ta(openOcean).^4 ...
    - par.sigma .* To(openOcean).^4;
end

