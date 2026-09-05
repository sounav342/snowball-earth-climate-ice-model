function par = integrate_EBM_1d_sphere(par,var2)
%% ice flow in 1 horizontal dimension using spherical coordinates.

set(0,'DefaultFigureVisible','off');

% Create Folder to Export Plots

Directory = string(pwd);
plotfolder_EBM = Directory+'/EBMFigures';

plotfolder = sprintf('%s/EBM-nonlinear-exp-%.2d-Q-%.2d-N-%.2d-eps-%.2d-CO2-%.2d-%s/',plotfolder_EBM,par.EBM_expnum,par.Qo,par.N,100*par.A_a,par.df_CO2*10,var2);
is_output_N = par.N==1 || floor(par.N/par.Nplot)*par.Nplot==par.N || par.N==par.Nt;
if is_output_N && ~exist(sprintf('%s',plotfolder),'file')
    [status,msg] = mkdir(sprintf('%s',plotfolder));
elseif is_output_N && exist(sprintf('%s',plotfolder),'file')
    [status,msg] = rmdir(sprintf('%s',plotfolder),'s');
    [status,msg] = mkdir(sprintf('%s',plotfolder));
end

% Create Folder to Export EBM Restarts (i.e. FIS output)
restartfolder_EBM = Directory+'/EBMRestart';
if ~exist(sprintf('%s',restartfolder_EBM),'file')
    [status,msg] = mkdir(sprintf('%s',restartfolder_EBM));
end

%% initialize variables:
Ta = zeros(par.nj,3);
To = zeros(par.nj,3);
qa = zeros(par.nj,3);
Ts = NaN(par.nj,3);
Sal = zeros(par.nj,3);
T_f = zeros(par.nj,3);
h = zeros(par.nj,3);

diff_to_tf=zeros(par.nj,3);

Sal_S_n = zeros(par.nj,1);

Qs_n = zeros(par.nj,1);

P_n = zeros(par.nj,1);
E_n = zeros(par.nj,1);

P_np1 = zeros(par.nj,1);
E_np1 = zeros(par.nj,1);

Sal_S_np1 = zeros(par.nj,1);

Qo_n = zeros(par.nj,1);
Io_n = zeros(par.nj,1);
Ho_n = zeros(par.nj,1);
G_n = zeros(par.nj,1);

Qo_np1 = zeros(par.nj,1);
Io_np1 = zeros(par.nj,1);
Ho_np1 = zeros(par.nj,1);
G_np1 = zeros(par.nj,1);

Is_n = zeros(par.nj,1);
Hs_n = zeros(par.nj,1);
Is_np1 = zeros(par.nj,1);
Hs_np1 = zeros(par.nj,1);

Qa_n = zeros(par.nj,1);
Ia_n = zeros(par.nj,1);
Ha_n = zeros(par.nj,1);

Qa_np1 = zeros(par.nj,1);
Ia_np1 = zeros(par.nj,1);
Ha_np1 = zeros(par.nj,1);

rhoa_n = zeros(par.nj,1);
rhoa_np1 = zeros(par.nj,1);

diff_To_n = zeros(par.nj,1);
diff_Ta_n = zeros(par.nj,1);
diff_qa_n = zeros(par.nj,1);
diff_sal_n = zeros(par.nj,1);

diff_To_np1 = zeros(par.nj,1);
diff_Ta_np1 = zeros(par.nj,1);
diff_qa_np1 = zeros(par.nj,1);
diff_sal_np1 = zeros(par.nj,1);

L_n = zeros(par.nj,1);
L_np1 = zeros(par.nj,1);

Fo_n = zeros(par.nj,1);
C_n = zeros(par.nj,1);
Pi_n = zeros(par.nj,1);

Fo_np1 = zeros(par.nj,1);
C_np1 = zeros(par.nj,1);
Pi_np1 = zeros(par.nj,1);

dhdt_frz = zeros(par.nj,1);
dhdt_melt = zeros(par.nj,1);

% Load Model Variables
[h,R,Ta,To,Ts,T_f,qa,Sal,h_o,iter_n,del_h_i,h_o_sal,Sal_init]=EBM_read_restart(par,var2);
fprintf(1,'running 1D spherical EBM @ N = %.2d, Q=%.2d, EBM_expnum=%.2d, epsilon =%.2d, CO2 =%.2d, ice %s\n',par.N,par.Qo,par.EBM_expnum,100*par.A_a,par.df_CO2*10,var2);

Ts(:,1) = fill_ice_surface_temperature(Ts(:,1),R(:,1),par);
Ts(:,2) = fill_ice_surface_temperature(Ts(:,2),R(:,2),par);

if par.N == 1

    Sal_init = Sal;

end

if (par.N == 1 && par.EBM_expnum == 22)
    df = 0;
else
    df = par.df_CO2;
end

par.h_initial = h(:,2);
par.R_initial = R(:,2);
par.To_initial = To(:,2);
par.To_init_ave = globmean(To(:,2),par);
par.Sal_initial = Sal(:,2);

interior_idx = (2:(par.nj-1)).';
s = par.s(:);
s_jphalf = 0.5 * (s(interior_idx) + s(interior_idx + 1));
s_jmhalf = 0.5 * (s(interior_idx) + s(interior_idx - 1));
inv_s = 1 ./ s(interior_idx);
dtheta_inv2 = 1 / par.dtheta_rad^2;
R2_inv = 1 / par.R^2;
qa_diff_scale = par.rho_a * par.h_q * par.D_q * R2_inv * dtheta_inv2;
Ta_diff_scale = par.rho_a * par.c_a * par.h_a * par.D_a * R2_inv * dtheta_inv2;
To_diff_scale = par.rho_o * par.c_o * h_o * R2_inv * dtheta_inv2;
Sal_diff_scale = h_o * par.D_o * R2_inv * dtheta_inv2;

dhdt_odiff_Sal = zeros(par.nj,1); % Sensible Heat is obv. zero when h=0

% Load n=1 Model Fields

Qs_n = Q_s(R(:,1),Ts(:,1),par);

P_n=P(Ta(:,1),Ts(:,1),qa(:,1),par,plotfolder,1);
E_n=E1(R(:,1),Ta(:,1),To(:,1),Ts(:,1),qa(:,1),par);

Qo_n = Q_o(R(:,1),par);
Io_n = I_o(R(:,1),Ta(:,1),To(:,1),df,par);
Ho_n = H_o(R(:,1),Ta(:,1),To(:,1),par);
G_n(:) = par.G;

Is_n = I_s(R(:,1),Ta(:,1),Ts(:,1),df,par);
Hs_n = H_s(R(:,1),Ta(:,1),Ts(:,1),par);

Qa_n = Q_a(R(:,1),par);
Ia_n = I_a(R(:,1),par,Ta(:,1),To(:,1),Ts(:,1),df,1);
Ha_n = H_a(R(:,1),Ho_n,Hs_n,par);

rhoa_n = rho_a(Ta(:,1),par);

C_n = Cflux(h(:,1),R(:,1),Ts(:,1),T_f(:,1),par);
Fo_n = F_o(C_n,G_n,par,R(:,1));
Pi_n = P_i(P_n,Ts(:,1),par);

Sal_S_n = Sal_source(dhdt_odiff_Sal,h(:,1),R(:,1),Ts(:,1),T_f(:,1),par);

L_n(:) = par.L_v;

track_melt = NaN(par.nt,1);
Toplot = NaN(par.nt+1,1);
T_fplot = NaN(par.nt+1,1);
Tsplot = NaN(par.nt+1,1);
Taplot = NaN(par.nt+1,1);
Eplot = NaN(par.nt+1,2);

Horz_D = linspace(par.D_o,0,90);

ngrid = numel(rhoa_np1);              % Total number of elements
grp = ceil((1:ngrid)'/45);       % Create group labels for each 45-element block
ice_line_idx = find(R(:,2) > 0);        % Find indices where the value is 1

if isempty(ice_line_idx)
    D_o = par.D_o;
else
    % For each group, find the maximum index (i.e. the last occurrence of 1)
    lastIndices = accumarray(grp(ice_line_idx), ice_line_idx, [max(grp), 1], @max, NaN);
    LI = lastIndices(1,1);
    if isnan(LI)
        D_o = par.D_o;
    else
        D_o = Horz_D(LI);
    end
end

T_f = (0.0901 - 0.0575 * Sal) - 7.61e-4 * (par.g * h * par.rho_i / 1e4) + 273.16;

Toplot(1) = globmean(To(:,2),par);
T_fplot(1) = globmean(T_f(:,2),par);
valid_Ts = ~isnan(Ts(:,2));
Tsplot(1) = sum(Ts(valid_Ts,2))/sum(valid_Ts);
Taplot(1) = globmean(Ta(:,2),par);
Eplot(1,:) = [0 0];

print = 1;

%% Time step equation:
for n=1:par.nt

    if n==print
        EBM_n = n;
        print = print + 18143;
    end

    % Update time step
    par.n=n;
    time_kyr=n*par.dt/par.year/1000;

    dhdt_odiff_Sal(interior_idx) = -par.Beta*(To(interior_idx,2) - T_f(interior_idx,2))/(par.rho_o*par.L_f); % + dhdt_geo;

    %% Advance qa, Ta, To, Ts, Sal in time using 2nd Order Adams Bashforth

    % Compute Model fields of current time step

    Qs_np1 = Q_s(R(:,2),Ts(:,2),par);

    P_np1=P(Ta(:,2),Ts(:,2),qa(:,2),par,plotfolder,n);

    E_np1=E1(R(:,2),Ta(:,2),To(:,2),Ts(:,2),qa(:,2),par);

    Qo_np1 = Q_o(R(:,2),par);
    Io_np1 = I_o(R(:,2),Ta(:,2),To(:,2),df,par);
    Ho_np1 = H_o(R(:,2),Ta(:,2),To(:,2),par);
    G_np1(:) = par.G;

    Is_np1 = I_s(R(:,2),Ta(:,2),Ts(:,2),df,par);
    Hs_np1 = H_s(R(:,2),Ta(:,2),Ts(:,2),par);

    Qa_np1 = Q_a(R(:,2),par);
    Ia_np1 = I_a(R(:,2),par,Ta(:,2),To(:,2),Ts(:,2),df,n);
    Ha_np1 = H_a(R(:,2),Ho_np1,Hs_np1,par);

    rhoa_np1 = rho_a(Ta(:,2),par);

    C_np1 = Cflux(h(:,2),R(:,2),Ts(:,2),T_f(:,2),par);
    Fo_np1 = F_o(C_np1,G_np1,par,R(:,2));
    Pi_np1 = P_i(P_np1,Ts(:,2),par);

    Sal_S_np1 = Sal_source(dhdt_odiff_Sal,h(:,2),R(:,2),Ts(:,2),T_f(:,2),par);

    L_np1(:) = par.L_v;

    if any(~isnan(Ts(:,3)))
        fprintf(1,'error: leading condition of Ts loop is wrong\n');
        error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
    end

    qa_lap_n = s_jphalf .* (qa(interior_idx + 1,1) - qa(interior_idx,1)) ...
        - s_jmhalf .* (qa(interior_idx,1) - qa(interior_idx - 1,1));
    qa_lap_np1 = s_jphalf .* (qa(interior_idx + 1,2) - qa(interior_idx,2)) ...
        - s_jmhalf .* (qa(interior_idx,2) - qa(interior_idx - 1,2));
    diff_qa_n(interior_idx) = qa_diff_scale .* inv_s .* qa_lap_n;
    diff_qa_np1(interior_idx) = qa_diff_scale .* inv_s .* qa_lap_np1;
    RHS_qa_n = (diff_qa_n(interior_idx) - (P_n(interior_idx) - E_n(interior_idx))) ./ (par.rho_a * par.h_q);
    RHS_qa_np1 = (diff_qa_np1(interior_idx) - (P_np1(interior_idx) - E_np1(interior_idx))) ./ (par.rho_a * par.h_q);
    qa(interior_idx,3) = qa(interior_idx,2) + par.dt * (1.5 * RHS_qa_np1 - 0.5 * RHS_qa_n);

    odiff_n = D_o .* (1 - R(interior_idx,1));
    odiff_np1 = D_o .* (1 - R(interior_idx,2));
    To_lap_n = s_jphalf .* (To(interior_idx + 1,1) - To(interior_idx,1)) ...
        - s_jmhalf .* (To(interior_idx,1) - To(interior_idx - 1,1));
    To_lap_np1 = s_jphalf .* (To(interior_idx + 1,2) - To(interior_idx,2)) ...
        - s_jmhalf .* (To(interior_idx,2) - To(interior_idx - 1,2));
    diff_To_n(interior_idx) = To_diff_scale .* odiff_n .* inv_s .* To_lap_n;
    diff_To_np1(interior_idx) = To_diff_scale .* odiff_np1 .* inv_s .* To_lap_np1;
    RHS_To_n = (diff_To_n(interior_idx) ...
        + (1 - R(interior_idx,1)) .* Qo_n(interior_idx) ...
        + (1 - R(interior_idx,1)) .* Io_n(interior_idx) ...
        + (1 - R(interior_idx,1)) .* Ho_n(interior_idx) ...
        + G_n(interior_idx) ...
        - par.L_v .* E_n(interior_idx) .* (1 - R(interior_idx,1)) ...
        - R(interior_idx,1) .* par.Beta .* (To(interior_idx,1) - T_f(interior_idx,1)) ...
        + (1 - R(interior_idx,1)) .* df) ./ (par.rho_o * par.c_o * h_o);
    RHS_To_np1 = (diff_To_np1(interior_idx) ...
        + (1 - R(interior_idx,2)) .* Qo_np1(interior_idx) ...
        + (1 - R(interior_idx,2)) .* Io_np1(interior_idx) ...
        + (1 - R(interior_idx,2)) .* Ho_np1(interior_idx) ...
        + G_np1(interior_idx) ...
        - par.L_v .* E_np1(interior_idx) .* (1 - R(interior_idx,2)) ...
        - R(interior_idx,2) .* par.Beta .* (To(interior_idx,2) - T_f(interior_idx,2)) ...
        + (1 - R(interior_idx,2)) .* df) ./ (par.rho_o * par.c_o * h_o);
    To(interior_idx,3) = To(interior_idx,2) + par.dt * (1.5 * RHS_To_np1 - 0.5 * RHS_To_n);

    Ta_lap_n = s_jphalf .* (Ta(interior_idx + 1,1) - Ta(interior_idx,1)) ...
        - s_jmhalf .* (Ta(interior_idx,1) - Ta(interior_idx - 1,1));
    Ta_lap_np1 = s_jphalf .* (Ta(interior_idx + 1,2) - Ta(interior_idx,2)) ...
        - s_jmhalf .* (Ta(interior_idx,2) - Ta(interior_idx - 1,2));
    diff_Ta_n(interior_idx) = Ta_diff_scale .* inv_s .* Ta_lap_n;
    diff_Ta_np1(interior_idx) = Ta_diff_scale .* inv_s .* Ta_lap_np1;
    RHS_Ta_n = (diff_Ta_n(interior_idx) + Qa_n(interior_idx) + Ia_n(interior_idx) + Ha_n(interior_idx) + par.L_v * P_n(interior_idx) + 2 * df) ...
        ./ (par.rho_a * par.c_a * par.h_a);
    RHS_Ta_np1 = (diff_Ta_np1(interior_idx) + Qa_np1(interior_idx) + Ia_np1(interior_idx) + Ha_np1(interior_idx) + par.L_v * P_np1(interior_idx) + 2 * df) ...
        ./ (par.rho_a * par.c_a * par.h_a);
    Ta(interior_idx,3) = Ta(interior_idx,2) + par.dt * (1.5 * RHS_Ta_np1 - 0.5 * RHS_Ta_n);

        %% Ocean Salinity, Sal

        sal_mean_n = globmean(Sal(:,1),par);
        sal_mean_np1 = globmean(Sal(:,2),par);

        Sal_lap_n = s_jphalf .* (Sal(interior_idx + 1,1) - Sal(interior_idx,1)) ...
            - s_jmhalf .* (Sal(interior_idx,1) - Sal(interior_idx - 1,1));
        Sal_lap_np1 = s_jphalf .* (Sal(interior_idx + 1,2) - Sal(interior_idx,2)) ...
            - s_jmhalf .* (Sal(interior_idx,2) - Sal(interior_idx - 1,2));
        diff_sal_n(interior_idx) = Sal_diff_scale .* inv_s .* Sal_lap_n;
        diff_sal_np1(interior_idx) = Sal_diff_scale .* inv_s .* Sal_lap_np1;
        RHS_sal_n = (diff_sal_n(interior_idx) ...
            + R(interior_idx,1) .* Sal_S_n(interior_idx) .* sal_mean_n ...
            - (1 - R(interior_idx,1)) .* ((P_n(interior_idx) - E_n(interior_idx)) / par.rho_o) .* sal_mean_n) ./ h_o;
        RHS_sal_np1 = (diff_sal_np1(interior_idx) ...
            + R(interior_idx,2) .* Sal_S_np1(interior_idx) .* sal_mean_np1 ...
            - (1 - R(interior_idx,2)) .* ((P_np1(interior_idx) - E_np1(interior_idx)) / par.rho_o) .* sal_mean_np1) ./ h_o;
        Sal(interior_idx,3) = Sal(interior_idx,2) + par.dt * (1.5 * RHS_sal_np1 - 0.5 * RHS_sal_n);

    RHS_Ts_n = (Qs_n(interior_idx) + Is_n(interior_idx) + Hs_n(interior_idx) - par.L_v * E_n(interior_idx) + C_n(interior_idx) + df) ...
        ./ (par.rho_i * par.c_i * par.d_i);
    RHS_Ts_np1 = (Qs_np1(interior_idx) + Is_np1(interior_idx) + Hs_np1(interior_idx) - par.L_v * E_np1(interior_idx) + C_np1(interior_idx) + df) ...
        ./ (par.rho_i * par.c_i * par.d_i);
    both_ice = R(interior_idx,1) > 0 & R(interior_idx,2) > 0;
    new_ice = R(interior_idx,1) == 0 & R(interior_idx,2) > 0;
    Ts(interior_idx,3) = NaN;
    Ts(interior_idx(both_ice),3) = Ts(interior_idx(both_ice),2) ...
        + par.dt * (1.5 * RHS_Ts_np1(both_ice) - 0.5 * RHS_Ts_n(both_ice));
    Ts(interior_idx(new_ice),3) = Ts(interior_idx(new_ice),2) + par.dt * RHS_Ts_np1(new_ice);

    h(interior_idx,3) = h(interior_idx,2);
    R(interior_idx,3) = R(interior_idx,2);

    T_f = (0.0901 - 0.0575 * Sal) - 7.61e-4 * (par.g * h * par.rho_i / 1e4) + 273.16;

    % set north and south boundary conditions
    Ta(1,3) = Ta(2,3);
    To(1,3) = To(2,3);
    qa(1,3) = qa(2,3);
    Ts(1,3) = Ts(2,3);
    Sal(1,3) = Sal(2,3);
    T_f(1,3) = T_f(2,3);

    Ta(par.nj,3) = Ta(par.nj-1,3);
    To(par.nj,3) = To(par.nj-1,3);
    qa(par.nj,3) = qa(par.nj-1,3);
    Ts(par.nj,3) = Ts(par.nj-1,3);
    Sal(par.nj,3) = Sal(par.nj-1,3);
    T_f(par.nj,3) = T_f(par.nj-1,3);

    h(1,3) = h(2,3);
    h(par.nj,3) = h(par.nj-1,3);
    R(1,3) = R(2,3);
    R(par.nj,3) = R(par.nj-1,3);

    if any(isnan(Ta(:,3))) || any(isnan(Ia_n)) || any(isnan(Ha_n))
        fprintf(1,'Check for NaNs\n');
        error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
    end

    % Check for Unphysical Values of Ts & Melt Surface Layer if Ts>Tf
    % melt water assumed to drain to ocean @ freezing temp

    track_melt(n) = globmean(dhdt_melt,par);

    % Check for Unphysical Ocean Temperatures
    % Once the EBM has equilibrated, freeze some ice if To<Tf or melt some ice (if it exists) where To>Tf

   % Check for Unphysical Ocean Temperatures
    if any(To(:,3)>373.15)
        fprintf(1,'error: To > Tboil\n');
        error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
    end

    if any(isnan(R(:,3))) || any(h(:,3)==10000)
        fprintf(1,'Check for NaNs in R or h = 1e4, i.e. resets\n');
        error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
    end

    [Eplot_n,netPrec,netEvap,Aplot,Surf,ToA] = ebalance(Ha_np1,par,R(:,2),Ta(:,2),To(:,2),Ts(:,2),T_f(:,2),dhdt_frz,dhdt_melt,n,C_np1,G_np1,P_np1,E_np1,var2,Qo_np1,Qs_np1,qa(:,2));
    Toplot(n+1) = globmean(To(:,2),par);
    T_fplot(n+1) = globmean(T_f(:,2),par);
    valid_Ts = ~isnan(Ts(:,2));
    Tsplot(n+1) = sum(Ts(valid_Ts,2))/sum(valid_Ts);
    Taplot(n+1) = globmean(Ta(:,2),par);
    Eplot(n+1,:) = Eplot_n;

    if n==par.nt
        par.Eplot_toa = Eplot(end,1);
        par.Eplot_sfc = Eplot(end,2);
        par.netPrec = netPrec;
        par.netEvap = netEvap;
        par.Aplot = Aplot;
    end

    st = ceil(par.nt/2);
    tt = (st:par.nt)/par.year;

    if is_output_N && n==par.nt

 end

    %% -----------
    %% plot Ta,To,qa,h,R,sal
    %% -----------
    if n==1 || floor(n/par.nplot)*par.nplot==n || n==par.nt || par.icelatlim==1 || par.icelatpole==1 || par.ebalance==1

        if is_output_N && (n==par.nt || par.icelatlim==1 || par.icelatpole==1 || par.ebalance==1)
            jb = 2; je = par.nj-1;

            prec_np1 = P_np1*par.year*100/par.rho_o;
            evap_np1 = E_np1*par.year*100/par.rho_o;

%             % Figure 1: Precipitation
%             % Figure 2: Evaporation
            % Figure 3: Net Solar Radiative Fluxes Absorbed by Surface
%             % Figure 4: Net Infrared Flux Absorbed by Surface
%             % Figure 5: Sensible Heat Fluxes Transferred to Surface

            % Figure 6: Surface Temperatures
            numyrs = (n/par.nt)*(par.Time/par.year); % time in yrs

            fig = figure('Visible','off');
            subplot(4,1,1)
            plot(90-par.theta(je:-1:jb),To(je:-1:jb,1)-par.T_f,'--ob',90-par.theta(je:-1:jb),To(je:-1:jb,2)-par.T_f,'--or',90-par.theta(je:-1:jb),To(je:-1:jb,3)-par.T_f,'--ok');

            yl = ylabel('C');
            xl = xlabel('Latitude');
            tl = title(sprintf('Ocean Temperature as a function of latitude at t=%.2d yrs',numyrs));
            set([gca,xl,yl,tl],'fontsize',10);
            subplot(4,1,2)
            plot(90-par.theta(je:-1:jb),Ts(je:-1:jb,1)-par.T_f,'--ob',90-par.theta(je:-1:jb),Ts(je:-1:jb,2)-par.T_f,'--or');

            yl = ylabel('C');
            xl = xlabel('Latitude');
            tl = title(sprintf('Ice Surface Temperature as a function of latitude at t=%.2d yrs',numyrs));
            set([gca,xl,yl,tl],'fontsize',10);

            subplot(4,1,3)
            plot(90-par.theta(je:-1:jb),Ta(je:-1:jb,1)-par.T_f,'--ob',90-par.theta(je:-1:jb),Ta(je:-1:jb,2)-par.T_f,'--or',90-par.theta(je:-1:jb),Ta(je:-1:jb,3)-par.T_f,'--ok');

            yl = ylabel('C');
            xl = xlabel('Latitude');
            tl = title(sprintf('Atmospheric Temperature as a function of latitude at t=%.2d yrs',numyrs));
            set([gca,xl,yl,tl],'fontsize',10);

            subplot(4,1,4)
            plot(90-par.theta(je:-1:jb),T_f(je:-1:jb,1)-par.T_f,'--ob',90-par.theta(je:-1:jb),T_f(je:-1:jb,2)-par.T_f,'--or',90-par.theta(je:-1:jb),T_f(je:-1:jb,3)-par.T_f,'--ok');
            yl = ylabel('C');
            xl = xlabel('Latitude');
            tl = title(sprintf('Freezing Temperature as a function of latitude at t=%.2d yrs',numyrs));
            set([gca,xl,yl,tl],'fontsize',10);

            if ishghandle(fig)
                try
                    saveas(fig,sprintf('%s/Tx-%.2d.png',plotfolder,n));
                    close(fig);
                catch
                    warning('integrate_EBM_1d_sphere:plotSave','Skipping Tx plot save at n=%d',n);
                end
            end

%             % Figure 7: Atmospheric Specific Humidity
%             % Figure 8: Fractional Field Height and Sea Glacier Thickness
%             % Figure 9: Diffusion

            % Figure 10: Ocean and Atmospheric Heat Transport
            % Total Ocean and Atmospheric Poleward Heat Transport in PW
            [atm_E,ocn_E,netatm_E,netocn_E,atm_diff,ocn_diff,netatm_diff,netocn_diff] = hflux(par,Surf,ToA,diff_To_np1,diff_Ta_np1,diff_qa_np1);

            % Figure 11: Ocean Salinity
            fig = figure('Visible','off');
            plot(90-par.theta(je:-1:jb),Sal(je:-1:jb,1),'--ob',90-par.theta(je:-1:jb),Sal(je:-1:jb,2),'--or',90-par.theta(je:-1:jb),Sal(je:-1:jb,3),'--ok');
            yl = ylabel('PPT');
            xl = xlabel('Latitude');
            tl = title('Oceanic Salinity');
            set([gca,xl,yl,tl],'fontsize',10);

            if ishghandle(fig)
                try
                    saveas(fig,sprintf('%s/Sal-%.2d.png',plotfolder,n));
                    close(fig);
                catch
                    warning('integrate_EBM_1d_sphere:plotSave','Skipping Sal plot save at n=%d',n);
                end
            end

            fig = figure('Visible','off');
            plot(90-par.theta(je:-1:jb),To(je:-1:jb,1)-par.T_f,'--ob',90-par.theta(je:-1:jb),T_f(je:-1:jb,2)-par.T_f,'--or');
            legend To Tf
            yl = ylabel('C');
            xl = xlabel('Latitude');
            tl = title('To-Tf');
            set([gca,xl,yl,tl],'fontsize',10);

            if ishghandle(fig)
                try
                    saveas(fig,sprintf('%s/To-Tf-%.2d.png',plotfolder,n));
                    close(fig);
                catch
                    warning('integrate_EBM_1d_sphere:plotSave','Skipping To-Tf plot save at n=%d',n);
                end
            end

        end

        Hcr = par.Hcr;

        S = (P_np1 - E_np1)/par.rho_i; % FIS requires units of m/s
        % Fo_np1 added to source term in the ice flow model
        P_E = (P_np1-E_np1)/par.rho_i;

        %% Write restart or export - XX
        if is_output_N && (n==par.nt || par.icelatlim==1 || par.icelatpole==1 || par.ebalance==1)

%             % Figure 10: Source Term
        end

        % prepare fields for export
        Ta(:,1) = Ta(:,2);
        Ta(:,2) = Ta(:,3);
        To(:,1) = To(:,2);
        To(:,2) = To(:,3);
        Ts(:,1) = Ts(:,2);
        Ts(:,2) = Ts(:,3);
        qa(:,1) = qa(:,2);
        qa(:,2) = qa(:,3);
        Sal(:,1) = Sal(:,2);
        Sal(:,2) = Sal(:,3);
        T_f(:,1) = T_f(:,2);
        T_f(:,2) = T_f(:,3);

        h(:,1)=h(:,2);
        h(:,2)=h(:,3);
        R(:,1)=R(:,2);
        R(:,2)=R(:,3);

        initial_cells = find(par.R_initial>0);
        if ~isempty(initial_cells)
            hx = zeros(numel(initial_cells),1);
            full_now = R(initial_cells,3)==1;
            partial_now = R(initial_cells,3)>0 & R(initial_cells,3)<1;
            hx(full_now) = par.h_initial(initial_cells(full_now));
            hx(partial_now) = par.R_initial(initial_cells(partial_now))*par.Hcr;
            par.h_init_ave = mean(hx);
        else
            par.h_init_ave = 0;
        end

        par.h_min = min(h(h(:,3)>0,3));
        par.h_max = max(h(:,3));

        par.h_end = h(:,3);
        par.R_end = R(:,3);

        % mean ocean temp
        par.To_ave = globmean(To(:,3),par);
        % max ocean temp
        par.To_max = max(To(:,3));
        % min ocean temp
        par.To_min = min(To(:,3));

        % Average h
        ice_cells = find(R(:,3)>0);
        if ~isempty(ice_cells)
            hx = h(ice_cells,3);
            partial_now = R(ice_cells,3)>0 & R(ice_cells,3)<1;
            hx(partial_now) = R(ice_cells(partial_now),3)*par.Hcr;
            par.h_ave = mean(hx);
        else
            par.h_ave = 0;
        end

% %             h([1:par.Nedge,par.Sedge:end],:) = globmean(balanced_h,par)

        v_sal = globmean(Sal(:,2),par).*h_o;
        Toplot = Toplot(1:n+1);
        T_fplot = T_fplot(1:n+1);
        Tsplot = Tsplot(1:n+1);
        Taplot = Taplot(1:n+1);
        Eplot = Eplot(1:n+1,:);
        track_melt = track_melt(1:n);
        save(sprintf('%s/restart-EBM-expnum-%.2d-Q-%.2d-resnum-%.2d-eps-%.2d-CO2-%.2d-%s.mat',restartfolder_EBM,par.EBM_expnum,par.Qo,par.N,100*par.A_a,par.df_CO2*10,var2),'h','R','Ta','To','Ts','T_f','Sal','h_o','qa','Hcr','S','P_E','diff_Ta_np1','diff_To_np1','diff_qa_np1','diff_to_tf','Toplot','P_np1','E_np1','T_fplot','dhdt_frz','dhdt_melt','v_sal','h_o_sal','Sal_init');

        if par.icelatlim==1 || par.icelatpole==1 || par.ebalance==1
            return;
        end
    end

    %% prepare for next time step:
    Ta(:,1) = Ta(:,2);
    Ta(:,2) = Ta(:,3);
    To(:,1) = To(:,2);
    To(:,2) = To(:,3);
    Ts(:,1) = Ts(:,2);
    Ts(:,2) = Ts(:,3);
    qa(:,1) = qa(:,2);
    qa(:,2) = qa(:,3);
    Sal(:,1) = Sal(:,2);
    Sal(:,2) = Sal(:,3);
    T_f(:,1) = T_f(:,2);
    T_f(:,2) = T_f(:,3);
    Qs_n = Qs_np1;

    P_n = P_np1;

    E_n = E_np1;

    Qo_n = Qo_np1;
    Io_n = Io_np1;
    Ho_n = Ho_np1;
    G_n = G_np1;

    Is_n = Is_np1;
    Hs_n = Hs_np1;

    Qa_n = Qa_np1;
    Ia_n = Ia_np1;
    Ha_n = Ha_np1;

    rhoa_n = rhoa_np1;

    L_n = L_np1;

    C_n = C_np1;
    Fo_n = Fo_np1;
    Pi_n = Pi_np1;
    Sal_S_n = Sal_S_np1;

    h(:,1)=h(:,2);
    h(:,2)=h(:,3);
    R(:,1)=R(:,2);
    R(:,2)=R(:,3);

    % Reset h,R,Ts,To
    h(:,3) = 10000;
    R(:,3) = NaN;
    Ts(:,3) = NaN;
    To(:,3) = 0;
    Ta(:,3) = 10000;
    qa(:,3) = 100;

    T_f(:,3) = 273.16;
    T_f(:,3) = par.T_f;

    % Reset Diagnostic Fields
    dhdt_frz = zeros(par.nj,1);
    dhdt_melt = zeros(par.nj,1);
end
