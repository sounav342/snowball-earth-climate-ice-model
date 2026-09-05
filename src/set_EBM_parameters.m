function par = set_EBM_parameters(EBM_expnum,N,Qo,var1,var2,CO2_forcing)
%SET_EBM_PARAMETERS Define physical, numerical, and scenario parameters.

par.alpha_o = 0.1; % ocean surface albedo (see A13.1)
par.alpha_i = 0.6; % check with Eli about ice albedo
par.alpha_a = 0.24;

par.Tau = 0.63; %Atmospheric transmissivity

par.df_CO2 = CO2_forcing;

par.N = N;
par.EBM_expnum = EBM_expnum;
par.Nt = 500;

par.Nplot = 1;

par.icelatlim = 0;
par.icelatpole = 0;
par.ebalance = 0;

par.A_a = var1;
par.version = sprintf('1d-sphere-partial-%d',EBM_expnum);

par.Beta = 400; %ocean-ice base heat flux coefficient, W/(m^2 K)
par.c_a = 1004.64; % specific heat of air, J/(kg K)
par.c_o = 4218; % specific heat of water, J/(kg K)
par.c_i = 2106; % specific heat of ice J/(kg K)

% C_D :  neutral drag coefficient (A13.3)
par.k = 0.4; % von Karman constant
par.v_a = 5; % ~40m level wind speed
par.z_a = 40; % reference height
par.z_s = 0.001; % surface roughness length for ice
par.z_o = 0.0001; % surface roughness length for ocean

par.CD_o = (par.k/(log(par.z_a/par.z_o)))^2;
par.CD_s = par.CD_o;

par.CD_evap = par.CD_o*1.65;
par.CD_sens = par.CD_o/4.5;

par.d_i = 0.5; % thickness of surf seasonal thermal layer (A13.4)
% Diffusivities (m^2 s^-1)
par.D_a = 1.3e6; % atmospheric heat diffusitity
par.D_q = 1.69e6; % atmospheric water vapor diffusivity  (A13.5)

par.D_o = 2.5e2; % ocean heat diffusivity

par.G = 0.06; % geothermal heat flux, W/m^2

par.h_a = 8400; % atmospheric thickness for heat, m
par.h_q = 6000; % atmospheric thickness for water vapor
par.h_o = 4000; % ocean mixed layer thickness

par.kappaI=1.36e-6; % heat cond. of ice

par.k_s = 0.2; % snow thermal conductivity, W/(mK)

% Latent Heat (J/kg)
par.L_f = 0.334e6; % latent heat of H20 fusion
par.L_v = 2.5104e6; % latent heat of H20 vaporization
par.L_s = 2.8440e6; % latent heat of H20 sublimation
% Densities (kg/m^3)
par.rho_a = 10^5/(287.04*287);
par.rho_o = 1000; % density of liquid water
par.rho_i = 900; % density of ice
par.rho_s = 250; % density of snow
par.sigma = 5.66961e-8; % Stefan-Boltzmann constant, W/(m^2 K^-4)
par.tau_p = 12*86400; % timescale for precipitation rate (s)
par.nn = 3; % Exponent of Glen's Law

par.k_i = 2.1; % ice thermal conductivity, W/(mK)

par.R_d = 287.058; % specific gas constant for dry air, J/(kg K)
par.R_v = 461.5; % specific gas constant for water vapor
par.epsilon = par.R_d/par.R_v;

par.Ta_o = 260; % 273
par.To_o = 290; % Average ocean temperature
par.Ts_o = 260; % Average surface temperature
par.p_a = 10^5; % Atmospheric pressure level
par.qa_o = 0.01; % Specific humidity

par.Sb = 50; %Initial Salinity in ppt

par.g=9.8; % gravity
par.T_f = 273.16; % melt freeze point of ice/water (K)

par.mu=par.rho_i/par.rho_o;
par.year=365*86400;
par.day=86400;

par.R=6300e3; % Earth Radius

% Thickness diffusivity, numerical only, make it as small as possible:
par.kappa=1;    % 1.0e0
par.kappa_2d=1.2e0;  % 1.0e0
par.T_surface_profile_type='warm';

par.tfac=20;

par.delta_T = 0.0001;

par.nplot=12*365*5; %48*7*54; % plot every x time steps
par.nplot_2d=100;

par.ni=89;
par.nj=89;
par.nk=42;
par.EQ = ceil(par.nj/2);

par.EBM_expnum=EBM_expnum;

scenarioFile = fullfile(fileparts(mfilename('fullpath')),'EBMInput', ...
    sprintf('exp_%02d.m',EBM_expnum));
if ~isfile(scenarioFile)
    error('Snowball:MissingScenario', ...
        'Experiment %02d is unavailable because %s is missing.',EBM_expnum,scenarioFile);
end
run(scenarioFile);

%% physical domain, including boundary points, is [2:nj-1,2:nk-1]
par.L = 2e7;
par.dzeta=1/(par.nk-3);
par.dx=par.L/(par.ni-3);
par.dy=par.L/(par.nj-3);

par.zeta =-par.mu+([1:par.nk]-2)*par.dzeta;
par.x=([1:par.ni]-2)*par.dx;
par.y=([1:par.nj]-2)*par.dy;

%% spherical coordinates:
par.theta_north=10;
par.dphi=360/(par.ni-3);
par.dtheta=(90-par.theta_north)*2/(par.nj-3);
par.dtheta1=(90-par.theta_north)*2/(par.ni-3);
par.dphi_rad=par.dphi*pi/180;
par.dtheta_rad=par.dtheta*pi/180;
par.phi=([1:par.ni]-2)*par.dphi;
par.theta=par.theta_north+([1:par.nj]-2)*par.dtheta;
par.theta1=par.theta_north+([1:par.ni]-2)*par.dtheta1;

par.s=sind(par.theta);
par.c=cosd(par.theta);
par.s_col = par.s(:);
par.c_col = par.c(:);
par.dc = diff(par.c_col);
par.dc_sum = sum(par.dc);
par.inv_R2_dtheta2 = 1/(par.R^2*par.dtheta_rad^2);
par.inv_2R_dtheta = 1/(2*par.R*par.dtheta_rad);

par.Time=par.tfac*par.year; % time to run the model for

par.dt = 1800;

par.nt=ceil(par.Time/par.dt); % par.nt = 2000

% Abbot et al., 2011
par.s2 = -0.482;
par.S = 1 - par.s2/2 + 1.5*par.s2*par.c.^2;
par.S_col = par.S(:);
par.Qo = Qo;

par.Q = par.Qo/globmean(par.S_col,par); % stellar constant normalized to the active forcing pattern

if par.N==1
    mkrestart(par,var2);
end

end
