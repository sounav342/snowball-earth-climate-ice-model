function Qs = Q_s(R,varargin)
% Qs: solar flux absorbed per unit ice-covered area.
% The grid-cell mean shortwave flux is R.*Qs + (1-R).*Qo.

if isempty(varargin)
    error('Q_s:missingParameters','Q_s requires R and par.');
end

par = varargin{end};
[~,alpha_i_local,~,R] = surface_albedo_fractional(R,par);

SW_down = par.Tau .* par.S_col(:) .* par.Q ./ 4;
Qs = zeros(par.nj,1);
iceCovered = R > 0;
Qs(iceCovered) = SW_down(iceCovered) .* (1 - alpha_i_local(iceCovered));
end
