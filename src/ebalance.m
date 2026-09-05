function [Eplot,netPrec,netEvap,Aplot,Surf,ToA] = ebalance( ...
    Ha_np1,par,R,Ta,To,Ts,T_f,dhdt_frz,dhdt_melt,n,C_np1,G_np1, ...
    P_np1,E_np1,~,Qo_np1,Qs_np1,~)
%EBALANCE Diagnose surface and top-of-atmosphere energy budgets.

R = validate_fractional_ice_cover(R,par,'ebalance');
Ta = Ta(:);
To = To(:);
Ts = Ts(:);
T_f = T_f(:);
Ha_np1 = Ha_np1(:);
dhdt_frz = dhdt_frz(:);
dhdt_melt = dhdt_melt(:);
C_np1 = C_np1(:);
G_np1 = G_np1(:);
P_np1 = P_np1(:);
E_np1 = E_np1(:);
Qo_np1 = Qo_np1(:);
Qs_np1 = Qs_np1(:);

Ts(isnan(Ts))=0;
% Net Shortwave (SW) at the surface
% Qs and Qo are per-area ice and ocean fluxes; area-weight them over R.
SW_net = (1 - R).*Qo_np1 + R.*Qs_np1;
% Net Longwave (LW) at the surface
LW_up = (1-R).*par.sigma.*To.^4 + R.*par.sigma.*Ts.^4;
LW_down = par.A_a*par.sigma*Ta.^4;
LW_net = LW_down - LW_up;
% Sensible Heat (SH) from the surface to the atmosphere
SH = Ha_np1;
% Latent Heat (LHe) from the surface to the atmosphere
LHe = E_np1*par.L_v;
% Conductive Flux (C) into the ice surface
C = (R).*C_np1;
% Geothermal Flux (G) into the ocean mixed layer (regardless of ice cover)
G = G_np1.*(1-R);
% Surface Melting or Basal Freezing (F)
F_melt = par.rho_i*par.L_f*dhdt_melt; % dhdt_melt<0
F_frz = par.rho_i*par.L_f*dhdt_frz;   % dhdt_frz>0
F_net = F_melt + F_frz;
% Sensible Heat from the ocean to the ice
% (zero in surface e bal for thick ice)
Rmask = R ~= 1;
SHoi = R.*par.Beta.*(To-T_f).*Rmask;
Surf = SW_net + LW_net - SH - LHe + C + G + F_net - SHoi;

% Net Shortwave (SWtoa) at the top of atmosphere
SWtoa = SW_net;
% Net Atmospheric Longwave (LWatm) to space
LWatm = par.A_a*par.sigma*Ta.^4;
% Net Surface Longwave (LWsfc) to space
LWsfc = (1-par.A_a)*LW_up;
% Latent Heat (LHp) of Condensation
% Residual P-E in the atmosphere

ToA = SWtoa - LWatm - LWsfc;

% Downwelling Solar Insolation (SW_down)
SW_down = par.S_col(:).*par.Q./4;
% Upwelling Solar Insolation (SW_up)
alpha_sfc = surface_albedo_fractional(R,par);
SW_up = alpha_sfc.*SW_down;
% Outgoing Longwave Radiation (OLR)
OLR = LWatm + LWsfc;

netSurf = globmean(Surf,par);
netToA = globmean(ToA,par);
netPrec = globmean(P_np1,par);
netEvap = globmean(E_np1,par);
netLW_up = globmean(LW_up,par);
netSW_down = globmean(SW_down,par);
netSW_up = globmean(SW_up,par);
netOLR = globmean(OLR,par);

Aplot = globmean(alpha_sfc,par);

Eplot = [netToA,netSurf];

if (par.N==1 || mod(par.N,par.Nplot)==0 || par.N==par.Nt) && n==par.nt
    fprintf(1,['Energy Balance: @n=par.nt, global mean albedo = %.2d, ' ...
        'netToA=%.2d & netSurf=%.2d & netPrec=%.2d & netEvap=%.2d & ' ...
        'netSW_down=%.2d & netSW_up=%.2d & netOLR=%.2d\n'], ...
        Aplot,netToA,netSurf,netPrec*100*par.year/par.rho_o, ...
        netEvap*100*par.year/par.rho_o,netSW_down,netSW_up,netOLR);
    if abs(netSurf)<0.5 && abs(netToA)<0.5
        fprintf(1,'Achieved Energy Balance\n');
    else
        fprintf(1,'Surface Energy Imbalance\n');
    end
end
