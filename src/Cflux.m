function C=Cflux(h,R,Ts,T_f,par)
% Conductive Flux underneath sea ice
C = zeros(par.nj,1);
Ts(isnan(Ts))=0;

full_ice = R == 1;
partial_ice = R ~= 0 & R ~= 1;

C(full_ice) = par.k_i .* (T_f(full_ice) - Ts(full_ice)) ./ h(full_ice);
C(partial_ice) = par.k_i .* (T_f(partial_ice) - Ts(partial_ice)) ./ par.Hcr;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Conductive Flux underneath sea ice
% % Ensure that To~=Tf under thick ice
% % if any(To ~= par.T_f & R == 1)
% % end
% % Find indices where R is 1
% % Find indices where R is not 1
% % Calculate conductive flux for R(j) == 1
% % Set conductive flux to 0 for R(j) == 0
% % Set conductive flux for R(j) == 2 using par.Hcr

