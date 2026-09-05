function par=set_FIS_parameters(EBM_expnum,N,Qo,CO2_forcing)

par.alpha_o = 0.1;
par.alpha_i = 0.6;
par.alpha_a = 0.24;

par.Tau = 0.63; %Atmospheric transmissivity
par.df_CO2 = CO2_forcing;

par.N = N;
par.EBM_expnum = EBM_expnum;
par.Nt = 500;

par.Nplot = 1;

par.Qo = Qo;

par.h_o = 4000; %xx 30

par.rho_i = 900;
par.c_i=2106;
par.kappaI=1.36e-6;

par.c_o = 4218;
par.G = 0.06; % geothermal heat flux, W/m^2
par.k_i = 2.1; % ice thermal conductivity, W/(m)

par.D_o = 2.5e2; % ocean heat diffusivity xx8e5

par.Beta = 400; %ocean-ice base heat flux coefficient, W/(m^2 K)

par.icelatlim = 0;
par.icelatpole = 0;
par.icelateq = 0;
par.icelatstable = 0;
par.icethicklim = 0;
par.noice = 0;

%% initial thickness if no restart file exists:
par.h0=1500;

par.g=9.8;
par.T_f=273.16;
par.rho_o=1000;
par.rho_i=900;
%% (J/kg, from wikipedia's latent_heat),
%% http://en.wikipedia.org/wiki/Latent_heat#Table_of_latent_heats :
par.L_f=334e3;
par.mu=par.rho_i/par.rho_o;
par.year=365*86400;
par.S0=0.015/par.year; % m/sec
par.L=2e7;
par.R=6300e3; % Earth radius
par.D=1000;
% Thickness diffusivity, numerical only, make it as small as possible:
par.kappa=1;    % 1.0e0
par.kappa_2d=1.2e0;  % 1.0e0
par.T_surface_profile_type='warm';

par.Time=1e6*par.year; % time to run the model for (n=2,51*par.year)
par.Time_2d=1e5*par.year;
par.dt=par.year;
par.dt_2d=0.5e1*par.year;
par.nt=ceil(par.Time/par.dt); % par.nt = 2000
par.nt_2d=ceil(par.Time_2d/par.dt_2d);

par.nplot=300;
par.nplot_2d=100;
par.do_use_imagesc=0;
par.plot_2d.min_h=NaN;
par.plot_2d.max_h=NaN;
par.plot_2d.do_h_ylabel=1;
par.nwrite_restart=min(par.nplot,floor(par.nt/10));

% Full-snowball FIS/EBM coupling threshold. Once ice is global, changes in
% thickness feed back through pressure melting, salinity, and conductive flux
% even when ice area no longer changes.
par.dh_ice_couple_abs = 10; % meters
par.dh_ice_couple_frac = 0.05;

par.ni=89;
par.nj=89;
par.nk=42;
par.EQ = ceil(par.nj/2);

% %% file with E-P/ melt-freeze forcing:

%% physical domain, including boundary points, is [2:nj-1,2:nk-1]
par.dzeta=1/(par.nk-3);
par.dx=par.L/(par.ni-3);
par.dy=par.L;
par.zeta =-par.mu+([1:par.nk]-2)*par.dzeta;
par.x=([1:par.ni]-2)*par.dx;
par.y=([1:par.nj]-2)*par.dy;

%% spherical coordinates:
par.theta_north=10;
par.dphi=360/(par.ni-3);
par.dtheta=(90-par.theta_north)*2/(par.nj-3);
par.dphi_rad=par.dphi*pi/180;
par.dtheta_rad=par.dtheta*pi/180;
par.phi=([1:par.ni]-2)*par.dphi;
par.theta=par.theta_north+([1:par.nj]-2)*par.dtheta;
par.s=sind(par.theta);
par.c=cosd(par.theta);
par.s_col = par.s(:);
par.c_col = par.c(:);
par.dc = diff(par.c_col);
par.dc_sum = sum(par.dc);
par.inv_R2_dtheta2 = 1/(par.R^2*par.dtheta_rad^2);
par.inv_2R_dtheta = 1/(2*par.R*par.dtheta_rad);
