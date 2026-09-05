function par = integrate_h_1d_sphere(par,var1,var2)
%% ice flow in 1 horizontal dimension using spherical coordinates.

set(0,'DefaultFigureVisible','off');

Directory = string(pwd);

% Create Folder to Export FIS Plots
plotfolder_FIS = Directory+'/FISFigures';
if ~exist(sprintf('%s',plotfolder_FIS),'file')
    [status,msg] = mkdir(sprintf('%s',plotfolder_FIS));
end

% Create Folder to Export FIS Restarts
restartfolder_FIS = Directory+'/FISRestart';
if ~exist(sprintf('%s',restartfolder_FIS),'file')
    [status,msg] = mkdir(sprintf('%s',restartfolder_FIS));
end

% Specify part of filename
fig_filename_1d='1d-sphere-nonlinear';

plotfolder = sprintf('%s/FIS-exp-%.2d-Q-%.2d-%s-N-%.2d-eps-%.2d-CO2-%.2d-%s/',plotfolder_FIS,par.EBM_expnum,par.Qo,fig_filename_1d,par.N,100*var1,par.df_CO2*10,var2);
is_output_N = par.N==1 || floor(par.N/par.Nplot)*par.Nplot==par.N || par.N==par.Nt;
plotfolder_run = sprintf('%s/FIS-exp-%.2d-Q-%.2d-%s-N-%.2d-eps-%.2d-CO2-%.2d-%s',plotfolder_FIS,par.EBM_expnum,par.Qo,fig_filename_1d,par.N,100*var1,par.df_CO2*10,var2);
if is_output_N && ~exist(plotfolder_run,'file')
    [status,msg] = mkdir(plotfolder_run);
elseif is_output_N && exist(plotfolder_run,'file')
    [status,msg] = rmdir(plotfolder_run,'s');
    [status,msg] = mkdir(plotfolder_run);
end

fprintf(1,'running FIS @ N = %.2d, EBM_expnum %.2d, Qo %.2d, eps %.2d, ice %s\n',par.N,par.EBM_expnum,par.Qo,100*var1,var2);

%% exponent of Glen's law.
nn=3;

%% initialize variables:
div_hv_n=NaN(par.nj,1);
div_hv_np1=NaN(par.nj,1);
kappa_del2_h_n=NaN(par.nj,1);
kappa_del2_h_np1=NaN(par.nj,1);
v_n=zeros(par.nj,1);
v_np1=zeros(par.nj,1); % vn plus 1

diff_sal_n = zeros(par.nj,1);

diff_sal_np1 = zeros(par.nj,1);

mask=ones(par.nj,1);
nan_mask=mask;
mask(1)=0; mask(par.nj)=0;
nan_mask(mask==0)=NaN;

% Load Model Variables

[h,R,Ta,To,Ts,T_f,qa,Hcr,S_init,Sal,h_o,h_o_sal,Sal_init]=FIS_read_restart(par,var1,var2);

par.Hcr = Hcr;

h_o_s = h_o_sal;

% Sensible Heat from Ocean to Ice Base in Source Function

index = 2:(par.nj-1);

dhdt_odiff = zeros(par.nj,1);
dhdt_odiff_Sal = zeros(par.nj,1);

dhdt_odiff(index) = -par.Beta*(To(index,2) - T_f(index,2))/(par.rho_i*par.L_f);% + dhdt_geo;

dhdt_odiff_Sal(index) = -par.Beta*(To(index,2) - T_f(index,2))/(par.rho_o*par.L_f);% + dhdt_geo;

% Assign initial Ts
T_surface = Ts(:,1:2);
T_ocean = To(:,1:2);

% Define Nedge_initial
arg_NH = find(h(1:par.EQ,2)==0,1,'first');
if ~isempty(arg_NH)
    if arg_NH==1
        par.Nedge_initial=NaN;
    elseif arg_NH>=2
        par.Nedge_initial=arg_NH-1;
    else
        fprintf(1,'error: edge\n');
        error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
    end
elseif sum(h(1:par.EQ,2)>0)==ceil(par.nj/2)
    % no h=0 found - ice covered
    par.Nedge_initial=0;
else
    fprintf(1,'unforseen Nedge routine\n');
    error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
end

% Define Sedge_initial
sh_len = length(h(par.EQ:end,2));
arg_SH = find(h(par.EQ:end,2)==0,1,'last');
if ~isempty(arg_SH)
    if arg_SH==sh_len
        par.Sedge_initial=NaN;
    elseif arg_SH<=(par.nj-1)
        par.Sedge_initial = arg_SH + par.EQ;
    else
        fprintf(1,'error: edge\n');
        error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
    end
elseif sum(h(par.EQ:end,2)>0)==ceil(par.nj/2)
    % no h=0 found - ice covered
    par.Sedge_initial = 0;
else
    fprintf(1,'unforseen Nedge routine\n');
    error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
end

par.h_initial = h(:,2);
par.R_initial = R(:,2);
par.To_initial = To(:,2);
par.Ta_initial = Ta(:,2);
par.To_init_ave = globmean(To(:,2),par);
par.Ta_init_ave = globmean(Ta(:,2),par);
par.Sal_init_ave = globmean(Sal(:,2),par);

% Checking Input File for Open Ocean
par = FIS_iceline(h,R,par,'start');

if par.noice==1 || par.icelatpole==1
    par.n = 0;
    par.h_end = h(:,2);
    par.R_end = R(:,2);
    par.h_max = max(h(:,2));
    positive_h = h(h(:,2)>0,2);
    if ~isempty(positive_h)
        par.h_min = min(positive_h);
    else
        par.h_min = 0;
    end
    ice_cells = find(R(:,2)>0);
    if ~isempty(ice_cells)
        hx = h(ice_cells,2);
        partial_now = R(ice_cells,2)>0 & R(ice_cells,2)<1;
        hx(partial_now) = R(ice_cells(partial_now),2)*par.Hcr;
        par.h_ave = mean(hx);
    else
        par.h_ave = 0;
    end
    par.To_ave = globmean(To(:,2),par);
    par.To_max = max(To(:,2));
    par.To_min = min(To(:,2));
    par.Ta_ave = globmean(Ta(:,2),par);
    par.Ta_max = max(Ta(:,2));
    par.Ta_min = min(Ta(:,2));
    par.Sal_ave = globmean(Sal(:,2),par);
    par.v_ave = 0;
    par.v_max = 0;
    par.v_min = 0;
    return;
end

print = 1;

h_init = h;

h_o_init = h_o;
s = par.s(:);
inv_R2_dtheta2 = 1/(par.R^2*par.dtheta_rad^2);
inv_2R_dtheta = 1/(2*par.R*par.dtheta_rad);

hhplot = NaN(par.nt+1,1);
hhplotmax = NaN(par.nt+1,1);
hhplotmin = NaN(par.nt+1,1);
hhplot(1) = globmean(h(:,2),par);
hhplotmax(1) = max(h(:,2));
hhplotmin(1) = min(h(:,2));

%% Time step h-equation:
for n=1:par.nt  %%%%%%%%%%%%%%%%%%%%%%%%%% for n = 1 or n = par.nt/2

    if n==print
        FIS_n = n;
        if n>1
            print = print + 100;
        else
            print = print + 99;
        end
    end

    % Update time step
    par.n=n;
    time_kyr=n*par.dt/par.year/1000;
    sgi = [];
    exsgi = [];
    isogci = [];
    oogci = [];
    sgilen = 0;
    exsgilen = 0;
    isogcilen = 0;
    oogcilen = 0;

    % Update v_n
    v_n = v_np1;

    % Update subgrid indices (sgi)
    % Locate (if any) new subgrid indices (existing subgrid cells in R)
    new_sg = h(:,2)<par.Hcr & h(:,2)>0;
    if any(new_sg)
        R(new_sg,2) = find_R(h(new_sg,2),par);
        h(new_sg,2) = 0;
    end

    % Update par.Nedge
    arg_NH = find(h(1:par.EQ,2)==0,1,'first');
    if ~isempty(arg_NH)
        if arg_NH==1
            par.Nedge=NaN;
        elseif arg_NH>=2
            par.Nedge=arg_NH-1;
        else
            fprintf(1,'error: edge\n');
            error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
        end
    elseif sum(h(1:par.EQ,2)>0)==ceil(par.nj/2)
        % no h=0 found - ice covered
        par.Nedge=0;
    else
        fprintf(1,'unforseen Nedge routine\n');
        error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
    end
    % Update par.Sedge
    sh_len = length(h(par.EQ:end,2));
    arg_SH = find(h(par.EQ:end,2)==0,1,'last');
    if ~isempty(arg_SH)
        if arg_SH==sh_len
            par.Sedge=NaN;
        elseif arg_SH<=(par.nj-1)
            par.Sedge = arg_SH + par.EQ;
        else
            fprintf(1,'error: edge\n');
            error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
        end
    elseif sum(h(par.EQ:end,2)>0)==ceil(par.nj/2)
        % no h=0 found - ice covered
        par.Sedge = 0;
    else
        fprintf(1,'unforseen Sedge routine\n');
        error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
    end

    % Check for invalid par.Nedge
    if (((par.Nedge<=2) || (par.Nedge>=par.EQ)) && (par.Nedge~=0))
        fprintf(1,'par.Nedge = %.2d\n',par.Nedge);
        fprintf(1,'*** par.Nedge invalid: advise stopping model\n');
        error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
    end
    % Check for invalid par.Sedge
    if (((par.Sedge<=par.EQ) || (par.Sedge>=(par.nj-1))) && (par.Sedge~=0))
        fprintf(1,'par.Sedge = %.2d\n',par.Sedge);
        fprintf(1,'*** par.Sedge invalid: advise stopping model\n');
        error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
    end

    % Identify subgrid indices (sgi) and extra subgrid indices (exsgi)
    all_sgi = find(R(:,2)>0 & R(:,2)<1);
    if ~isempty(all_sgi)
        if par.Nedge==0 && par.Sedge==0
            sgilen=0;
            exsgi = all_sgi;
            exsgilen = length(exsgi);
            fprintf(1,'this an error?\n');
            error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
        elseif par.Nedge==0
            sgi = all_sgi(all_sgi==(par.Sedge-1));
            sgilen = length(sgi);
            if sgilen==0
                exsgi = all_sgi;
                exsgilen = length(exsgi);
            else
                exsgi = all_sgi(1:(end-1));
                exsgilen = length(exsgi);
            end
        elseif par.Sedge==0
            sgi = all_sgi(all_sgi==(par.Nedge+1));
            sgilen = length(sgi);
            if sgilen==0
                exsgi = all_sgi;
                exsgilen = length(exsgi);
            else
                exsgi = all_sgi(2:end);
                exsgilen = length(exsgi);
            end
        else
            sgi = all_sgi(all_sgi==(par.Nedge+1) | all_sgi==(par.Sedge-1));
            sgilen = length(sgi);

            if sgilen==0
                exsgi = all_sgi;
                exsgilen = length(exsgi);
            elseif sgilen==1
                if sgi>par.Nedge && sgi<=par.EQ
                    exsgi = all_sgi(2:end);
                    exsgilen = length(exsgi);
                elseif sgi<par.Sedge && sgi>=par.EQ
                    exsgi = all_sgi(1:(end-1));
                    exsgilen = length(exsgi);
                else
                    fprintf(1,'error: sgi calc - two hemis\n');
                    error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
                end
            elseif sgilen==2
                exsgi = all_sgi(2:(end-1));
                exsgilen = length(exsgi);
            end
        end
        v_np1(:) = 0; % XX
    end

    %% Check new sgi
    if sgilen==2
        if all(sgi < par.EQ)
            fprintf(2,'error: multiple SCs (not eSCs) in one hemisphere');
            error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
        end
        if all(sgi > par.EQ)
            fprintf(2,'error: multiple SCs (not eSCs) in one hemisphere');
            error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
        end
    end

    % Check for invalid NH subgrid-terminus regime
    if ~isempty(all_sgi) && par.Nedge==0 && min(all_sgi)<=par.EQ
        fprintf(1,'par.Nedge = %.2d; min(sgi)=%.2d\n',par.Nedge,min(all_sgi));
        fprintf(1,'*** invalid NH subgrid-terminus regime: advise stopping model\n');
    end
    % Check for invalid SH subgrid-terminus regime
    if ~isempty(all_sgi) && par.Sedge==0 && max(all_sgi)>=par.EQ
        fprintf(1,'par.Sedge = %.2d; max(sgi)=%.2d\n',par.Sedge,max(all_sgi));
        fprintf(1,'*** invalid SH subgrid-terminus regime: advise stopping model\n');
    end

    %% Edited isogci and oogci to work when there aren't two terminuses
    % Identify isolated grid cells (isogci)
    all_gci = find(R(:,2)==1);
    if par.Nedge==0 && par.Sedge~=0
        minval = par.EQ+1;
        maxval = par.Sedge-1;
        if minval > maxval
            isogcilen = 0;
        else
            indexrange = all_gci > minval & all_gci < maxval;
            isogci = all_gci(indexrange);
            isogcilen = length(isogci);
        end
    elseif par.Sedge==0 && par.Nedge~=0
        minval = par.Nedge+1;
        maxval = par.EQ-1;
        if minval > maxval
            isogcilen = 0;
        else
            indexrange = all_gci > minval & all_gci < maxval;
            isogci = all_gci(indexrange);
            isogcilen = length(isogci);
        end
    else
        minval = par.Nedge+1;
        maxval = par.Sedge-1;
        indexrange = all_gci > minval & all_gci < maxval;
        isogci = all_gci(indexrange);
        isogcilen = length(isogci);
    end

    % Identify Open Ocean Grid Cells (oogci)
    all_ogi = find(R(:,2)==0);
    if ~isempty(all_ogi)
        if par.Nedge==0 && par.Sedge~=0
            minval = par.EQ+1;
            maxval = par.Sedge-1;
            if minval > maxval
                oogcilen = 0;
            else
                indexrange = all_ogi > minval & all_ogi < maxval;
                oogci = all_ogi(indexrange);
                oogcilen = length(oogci);
            end
        elseif par.Nedge~=0 && par.Sedge==0
            minval = par.Nedge+1;
            maxval = par.EQ-1;
            if minval > maxval
                oogcilen = 0;
            else
                indexrange = all_ogi > minval & all_ogi < maxval;
                oogci = all_ogi(indexrange);
                oogcilen = length(oogci);
            end
        else
            minval = par.Nedge+1;
            maxval = par.Sedge-1;
            indexrange = all_ogi > minval & all_ogi < maxval;
            oogci = all_ogi(indexrange);
            oogcilen = length(oogci);
        end
    else
        oogcilen = 0;
    end

    % Update the operable domain
    if ((par.Nedge==0)&(par.Sedge==0))
        par.domain = [2:par.nj-1];
    elseif ((par.Nedge~=0)&(par.Sedge==0))
        par.domain = [2:par.Nedge,ceil(par.nj/2):par.nj-1];
    elseif ((par.Nedge==0)&(par.Sedge~=0))
        par.domain = [2:ceil(par.nj/2),par.Sedge:par.nj-1];
    elseif ((par.Nedge~=0)&(par.Sedge~=0))
        par.domain = [2:par.Nedge,par.Sedge:par.nj-1];
    end

    % Extrapolate Surface Temperature of the Ice, if needed
    T_surface(:,2) = Ts_extrap_FIS(h(:,2),R,T_surface(:,2),par);

%     Extrapolate Mixed Layer Ocean Temperature, if needed

    % Transition from Partial to Global requires default v_np1 field
    if (par.Nedge==0 && par.Sedge==0) && any(isnan(v_n))
        v_np1(:) = 0;
    end

    % Update Source Function
  [S,dhdt_cond] = update_S(dhdt_odiff,h(:,2),R(:,2),S_init,T_surface(:,2),T_f(:,2),par);

    if n==1;
        % this requires same h,R initial conditions for n=1:2
        % in other words, the initial conditions should be equilibrium fields
        v_n=v_np1;
        B = ones(par.nj,1);
    end

    % Edit v-np1 & B if advance of ice produces new grid cell for V Loop
    % This step is important because initialization v-np1 is the
    % same as v-n, which has a smaller real-valued domain.
    % This step is only an issue if v_def or v_undef are <3 or >par.nj-2

    % New subgrids have a B, but no v
    v_undef = find(R(:,2)>0 & R(:,2)<1 & ...
        ~isnan(B) & ~isnan(v_np1));
    if ~isempty(v_undef)
        minval = 3;
        maxval = par.nj-2;
        indexrange = v_undef>=minval & v_undef<=maxval;
        v_undef = v_undef(indexrange);
        v_np1(v_undef) = NaN; % uninitalize v in new sg
        B(v_undef) = 1; % initialize B in new sg
    end

    % New grids have a B and a v
    v_def = find(h(:,2)>0 & ...
        (isnan(B) | isnan(v_np1)));
    if ~isempty(v_def)
        minval = 3;
        maxval = par.nj-2;
        indexrange = v_def>=minval & v_def<=maxval;
        v_def = v_def(indexrange);
        v_np1(v_def) = 0; % initalize v in new grid
        B(v_def) = 1; % initalize B in new grid
    end

    % Equating # of NaNs in v_n and v_np1 for Height Loop
    % Note: If v_n has more NaNs than v_np1, then the height loop below will
    % iterate over NaNs in v_n while iterating over real number values in v_np1.
    if sum(isnan(v_n))>sum(isnan(v_np1))
        fill_idx = isnan(v_n(3:(par.nj-2))) & ~isnan(v_np1(3:(par.nj-2)));
        fill_pos = find(fill_idx) + 2;
        v_n(fill_pos) = v_np1(fill_pos);
    end

    %% calculate velocity:
    %% -------------------
    max_iter_eff_viscosity=1000;
    par.iter_eff_viscosity=0;
    diff_eff_viscosity=1;
    while par.iter_eff_viscosity<=max_iter_eff_viscosity ...
            && diff_eff_viscosity>1.e-4;
        %% Average viscosity over relevant temperature range:
        B_old=B;

        [B,B_tN,B_tS,B_tSd_ftheta_v_N,d_ftheta_v_S]=calc_eff_viscosity_1d_sphere(n,par,T_f(:,2),T_surface(:,2),h(:,2),v_np1,nn);

        % Flow vs No-Flow Switch
        if strcmp(var2,'flow')==1
            [v_np1,par]=solve_v_1d_matrix_form_sphere(n,par,B,B_tN,B_tS,squeeze(h(:,2)));
        elseif strcmp(var2,'no-flow')==1
            v_np1 = zeros(par.nj,1);
            vmask = B; vmask(~isnan(vmask))=1;
            v_np1 = v_np1.*vmask;
        else
            fprintf(1,'error: flow condition unspecified\n');
            error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
        end

        max_B=max(B(mask~=0));
        diff_eff_viscosity=max(abs(B-B_old))/max_B;
        par.iter_eff_viscosity=par.iter_eff_viscosity+1;
        if par.iter_eff_viscosity==max_iter_eff_viscosity
            fprintf(1,'*** at n=%d, reached max_iter_eff_viscosity=%d at n=%d, with diff_eff_viscosity=%g\n' ...
                ,n,max_iter_eff_viscosity,par.n,diff_eff_viscosity);
        end
    end

    active_v_idx = par.domain(:);
    if any(isnan(v_np1(active_v_idx)))
        fprintf(2,'***error: NaN velocity inside active ice domain\n');
        error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
    end

    %% advance h in time using second order Adams Bashforth:
    %% NOTE: code currently assumes h does not go below a certain minimum

    %% Reinitialize variables:
    div_hv_n=NaN(par.nj,1);
    div_hv_np1=NaN(par.nj,1);
    kappa_del2_h_n=NaN(par.nj,1);
    kappa_del2_h_np1=NaN(par.nj,1);
    RHS_n = NaN(par.nj,1);
    RHS_np1 = NaN(par.nj,1);
    interior_domain = par.domain(par.domain ~= par.Nedge & par.domain ~= par.Sedge);
    if ~isempty(interior_domain)
        s_jmhalf = 0.5 * (s(interior_domain) + s(interior_domain - 1));
        s_jphalf = 0.5 * (s(interior_domain) + s(interior_domain + 1));

        kappa_del2_h_n(interior_domain) = par.kappa .* ...
            ((s_jphalf .* (h(interior_domain + 1,1) - h(interior_domain,1)) .* mask(interior_domain) .* mask(interior_domain + 1) ...
            - s_jmhalf .* (h(interior_domain,1) - h(interior_domain - 1,1)) .* mask(interior_domain) .* mask(interior_domain - 1)) ...
            .* inv_R2_dtheta2 ./ s(interior_domain));
        div_hv_n(interior_domain) = ...
            (s(interior_domain + 1) .* h(interior_domain + 1,1) .* v_n(interior_domain + 1) ...
            - s(interior_domain - 1) .* h(interior_domain - 1,1) .* v_n(interior_domain - 1)) ...
            .* inv_2R_dtheta ./ s(interior_domain);
        RHS_n(interior_domain) = -div_hv_n(interior_domain) + kappa_del2_h_n(interior_domain) + S(interior_domain);

        kappa_del2_h_np1(interior_domain) = par.kappa .* ...
            ((s_jphalf .* (h(interior_domain + 1,2) - h(interior_domain,2)) .* mask(interior_domain) .* mask(interior_domain + 1) ...
            - s_jmhalf .* (h(interior_domain,2) - h(interior_domain - 1,2)) .* mask(interior_domain) .* mask(interior_domain - 1)) ...
            .* inv_R2_dtheta2 ./ s(interior_domain));
        div_hv_np1(interior_domain) = ...
            (s(interior_domain + 1) .* h(interior_domain + 1,2) .* v_np1(interior_domain + 1) ...
            - s(interior_domain - 1) .* h(interior_domain - 1,2) .* v_np1(interior_domain - 1)) ...
            .* inv_2R_dtheta ./ s(interior_domain);
        RHS_np1(interior_domain) = -div_hv_np1(interior_domain) + kappa_del2_h_np1(interior_domain) + S(interior_domain);
        h(interior_domain,3) = h(interior_domain,2) + par.dt .* (1.5 .* RHS_np1(interior_domain) - 0.5 .* RHS_n(interior_domain));
    end

    edge_domain = par.domain(par.domain == par.Nedge | par.domain == par.Sedge);
    for i = 1:length(edge_domain)
        j = edge_domain(i);

        if j == par.Nedge
            s_jmhalf=0.5*(s(j)+s(j-1));
            s_jphalf=0.5*(s(j)+s(j+1));
            kappa_del2_h_n(j)=par.kappa*( ...
                (1/(par.R^2*s(j)))* ...
                (-s_jmhalf*(h(j,1)-h(j-1,1))*mask(j-1)*mask(j-2) ...
                )/par.dtheta_rad^2);

            div_hv_n(j)= (1/(par.R*s(j)))*(...
                (s(j)*h(j,1)*v_n(j)-s(j-1)*h(j-1,1)*v_n(j-1)) ...
                /(par.dtheta_rad) ...
                );
            RHS_n=-div_hv_n(j)+kappa_del2_h_n(j)+S(j);

            kappa_del2_h_np1(j)=par.kappa*( ...
                (1/(par.R^2*s(j)))* ...
                (-s_jmhalf*(h(j,2)-h(j-1,2))*mask(j-1)*mask(j-2) ...
                )/par.dtheta_rad^2);

            div_hv_np1(j)= (1/(par.R*s(j)))*(...
                (s(j)*h(j,2)*v_np1(j)-s(j-1)*h(j-1,2)*v_np1(j-1)) ...
                /(par.dtheta_rad) ...
                );
            RHS_np1=-div_hv_np1(j)+kappa_del2_h_np1(j)+S(j);

            h(j,3)=h(j,2)+par.dt*(1.5*RHS_np1-0.5*RHS_n);

        elseif j == par.Sedge

            s_jmhalf=0.5*(s(j)+s(j-1));
            s_jphalf=0.5*(s(j)+s(j+1));
            kappa_del2_h_n(j)=par.kappa*( ...
                (1/(par.R^2*s(j)))* ...
                (s_jphalf*(h(j+1,1)-h(j,1))*mask(j+2)*mask(j+1) ...
                )/par.dtheta_rad^2);

            div_hv_n(j)= (1/(par.R*s(j)))*(...
                (s(j+1)*h(j+1,1)*v_n(j+1)-s(j)*h(j,1)*v_n(j)) ...
                /(par.dtheta_rad) ...
                );
            RHS_n=-div_hv_n(j)+kappa_del2_h_n(j)+S(j);

            kappa_del2_h_np1(j)=par.kappa*( ...
                (1/(par.R^2*s(j)))* ...
                (s_jphalf*(h(j+1,2)-h(j,2))*mask(j+2)*mask(j+1) ...
                )/par.dtheta_rad^2);

            div_hv_np1(j)= (1/(par.R*s(j)))*(...
                (s(j+1)*h(j+1,2)*v_np1(j+1)-s(j)*h(j,2)*v_np1(j)) ...
                /(par.dtheta_rad) ...
                );
            RHS_np1=-div_hv_np1(j)+kappa_del2_h_np1(j)+S(j);

            h(j,3)=h(j,2)+par.dt*(1.5*RHS_np1-0.5*RHS_n);

        end

        % Final Assignment of h and R fields
        if h(j,3)<=0
            % Ice Melts to Open Ocean
            h(j,3) = 0;
            R(j,3) = 0;
            % xx - additional heat input to ocean?
        elseif h(j,3)<par.Hcr && h(j,3)>0
            R(j,3) = find_R(h(j,3),par);
            h(j,3) = 0;
        else
            % grid scale ice thickness change
            R(j,3) = find_R(h(j,3),par);
        end

    end

    active_domain = par.domain(:);
    melt_idx = active_domain(h(active_domain,3) <= 0);
    subgrid_idx = active_domain(h(active_domain,3) > 0 & h(active_domain,3) < par.Hcr);
    full_idx = active_domain(h(active_domain,3) >= par.Hcr);
    h(melt_idx,3) = 0;
    R(melt_idx,3) = 0;
    R(subgrid_idx,3) = find_R(h(subgrid_idx,3),par);
    h(subgrid_idx,3) = 0;
    R(full_idx,3) = find_R(h(full_idx,3),par);

    %%% Grid and Subgrid Parametrization: Advance and Retreat %%%
    if (sgilen==0 && (par.Nedge==0 && par.Sedge==0))
        % global ice cover - no advection at ice edge
        % contintuity eq solved in applicable domain
        % WARNING: if you update R here then you erase new sgs!
    elseif ((sgilen==0) && ((par.Nedge~=0) || (par.Sedge~=0)))
        % partial ice cover, no subgrid
        if ((par.Nedge==(par.EQ-1)) && (par.Sedge==(par.EQ+1)))
            % advection into EQ ocean cell on grid scale
            ogi = par.EQ;
            [h,R,par] = ogrid(ogi,h,R,v_np1,par,S);
        elseif par.Nedge==0
            % NH ice shelf reaches EQ, advection in SH
            ogi = par.Sedge-1;
            [h,R,par] = ogrid(ogi,h,R,v_np1,par,S);
        elseif par.Sedge==0
            % SH ice shelf reaches EQ, advection in NH
            ogi = par.Nedge+1;
            [h,R,par] = ogrid(ogi,h,R,v_np1,par,S);
        else
            % advection in NH & SH
            ogi = [par.Nedge+1,par.Sedge-1];
            [h,R,par] = ogrid(ogi,h,R,v_np1,par,S);
        end
    elseif sgilen==1
        % partial ice cover, 1 subgrid cell
        if par.Nedge==0
            % NH ice shelf reaches EQ, subgrid advection in SH
            % check sgi = par.Sedge-1
            [h,R,par] = sgrid(sgi,h,R,v_np1,par,S,n);
        elseif par.Sedge==0
            % SH ice shelf reaches EQ, subgrid in NH
            % check sgi = par.Nedge+1
            [h,R,par] = sgrid(sgi,h,R,v_np1,par,S,n);
        else
            if sgi < par.EQ
                % subgrid advection in NH, grid advection in SH
                ogi = par.Sedge-1;
                [h,R,par] = ogrid(ogi,h,R,v_np1,par,S);
                [h,R,par] = sgrid(sgi,h,R,v_np1,par,S,n);
            elseif sgi > par.EQ
                % subgrid advection in SH, grid advection in NH
                ogi = par.Nedge+1;
                [h,R,par] = ogrid(ogi,h,R,v_np1,par,S);
                [h,R,par] = sgrid(sgi,h,R,v_np1,par,S,n);
            elseif sgi == par.EQ
                % subgrid advection into EQ from NH and SH
                % check N = EQ-1, S = EQ+1
                [h,R,par] = sgrid(sgi,h,R,v_np1,par,S,n);
            end
        end
    elseif sgilen==2
        % partial ice cover, 2 subgrid cells
        [h,R,par] = sgrid(sgi,h,R,v_np1,par,S,n);
    end

    %% User's Note:  Mass Balance of Extra Subgrids is determined
    %% by the source function. If a hemisphere has more than two subgrids, the
    %% subgrids more than two or more cells away from the terminus are referred to
    %% as 'extra subgrids'. In the first step, the continuity equation is
    %% solved in the subgrid cells that border each terminus. If the terminus isn't
    %% adjacent to a subgrid, then the ogrid function calculates the advection of
    %% ice into the adjacent open ocean cell.

    %%% Mass Balance of Extra Subgrids %%%
    if exsgilen>0
        [h,R] = exsgrid(h,par,R,S,exsgi,v_np1,n);
    end

    %%% Mass Balance of Isolated Ice Shelves (Grid Cells) %%%
    if isogcilen>0
        [h,R] = isogrid(h,par,R,S,isogci,n);
    end

    %%% Mass Balance of Open Ocean Grid Cell %%%
    % Different class of ogi, which is the adjacent ocean cell to a terminus.
    % Open Ocean Grid Cell (oogci) defined as an ocean cell not adjacent to a terminus.
    if oogcilen>0
        [h,R] = OOgrid(h,par,R,S,oogci);
    end

    %% set north and south boundary conditions of h,R:
    h(1,3)=h(2,3);
    h(par.nj,3)=h(par.nj-1,3);
    R(1,3) = R(2,3);
    R(par.nj,3) = R(par.nj-1,3);

    % Check for NaNs in h,R
    if any(isnan(h(:,3))) || any(isnan(R(:,3)))
        fprintf(1,'checkpoint why NaN\n');
        error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
    end

    %% check for non-physical values of h & R:
    found_negative_h=0;
    if min(h(:,3))<0 || min(R(:,3))<0 || max(R(:,3))>1
        fprintf(1,'min(h)=%g; max(h)=%g; min(R)=%g; max(R)=%g;',...
            min(h(:,3)),max(h(:,3)),min(R(:,3)),max(R(:,3)));
        fprintf(1,'*** stopping due to h<0 or R<0, n=%d ***\n',par.n);
        found_negative_h=1;
    end

%     Check for Unphysical Values of To and Freeze Ice if To<Tf
%      or melt ice if To>Tf under Thick Ice
    % Update par for Equilibrium Ice Line Checks
    par = FIS_iceline(h,R,par,'end');

    is_fis_stop = found_negative_h==1 || par.icelatlim==1 || par.icelatpole==1 ...
        || par.icelateq==1 || par.icethicklim==1 || par.icelatstable==1;
    hhplot(n+1) = globmean(h(:,2),par);
    hhplotmax(n+1) = max(h(:,2));
    hhplotmin(n+1) = min(h(:,2));

    %% -----------
    %% plot h,u,v:
    %% -----------
    if is_output_N && (n==1 || n==par.nt || is_fis_stop)
        write_FIS_diagnostic_plots(plotfolder,n,par,h(:,2),R(:,2),T_surface(:,2),To(:,2), ...
            S,dhdt_cond,dhdt_odiff,S_init,div_hv_n,kappa_del2_h_n,v_np1,nan_mask);

        if is_fis_stop && n<par.nt && found_negative_h==0
            T_surface_end = T_surface(:,2);
            T_surface_end(R(:,3)==0) = NaN;
            T_surface_end = fill_ice_surface_temperature(T_surface_end,R(:,3),par);
            [S_end,dhdt_cond_end] = update_S(dhdt_odiff,h(:,3),R(:,3),S_init,T_surface_end,T_f(:,2),par);
            write_FIS_diagnostic_plots(plotfolder,n+1,par,h(:,3),R(:,3),T_surface_end,To(:,2), ...
                S_end,dhdt_cond_end,dhdt_odiff,S_init,div_hv_np1,kappa_del2_h_np1,v_np1,nan_mask);
        end

        if found_negative_h==1
            return
        end
    end %% plotting

    %% Write restart
    % Caution: this code rewrites over existing fields for export
    if n==par.nt || par.icelatlim==1 || par.icelatpole==1 || par.icelateq==1 || par.icethicklim==1 || par.icelatstable==1

        theta=par.theta;
        phi=par.phi;

        % Update h(:,3),R(:,3) before writing Restart
        new_sg = h(:,3)<par.Hcr & h(:,3)>0;
        if any(new_sg)
            R(new_sg,3) = find_R(h(new_sg,3),par);
            h(new_sg,3) = 0;
        end

        % Output FISEBM time-series fields
        par.h_min = min(h(h(:,3)>0,3));
        par.h_max = max(h(:,3));

        ice_cells = find(R(:,3)>0);
        if ~isempty(ice_cells)
            hx = h(ice_cells,3);
            partial_now = R(ice_cells,3)>0 & R(ice_cells,3)<1;
            hx(partial_now) = R(ice_cells(partial_now),3)*par.Hcr;
            par.h_ave = mean(hx);
        else
            par.h_ave = 0;
        end

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

        % latitudinal ice extent
        par.h_end = h(:,3);
        par.R_end = R(:,3);

        % mean ocean temp
        par.To_ave = globmean(To(:,3),par);
        % max ocean temp
        par.To_max = max(To(:,3));
        % min ocean temp
        par.To_min = min(To(:,3));

        % mean atmos temp
        par.Ta_ave = globmean(Ta(:,3),par);
        % max atmos temp
        par.Ta_max = max(Ta(:,3));
        % min atmos temp
        par.Ta_min = min(Ta(:,3));

        %mean salinity
        par.Sal_ave = globmean(Sal(:,3),par);

        % mean velocity
        valid_v = ~isnan(v_np1);
        if any(valid_v)
            par.v_ave = mean(abs(v_np1(valid_v)))*par.year;
        else
            par.v_ave = 0;
        end

        par.v_max = max(abs(v_np1))*par.year;
        par.v_min = min((v_np1(v_np1>0)))*par.year;

        par.v_min1 = min(abs(v_np1(v_np1>0)))*par.year;

        h(:,2) = h(:,3);
        R(:,2) = R(:,3);
        T_f(:,2) = T_f(:,3);

        T_surface(:,1) = fill_ice_surface_temperature(T_surface(:,1),R(:,1),par);
        T_surface(:,2) = fill_ice_surface_temperature(T_surface(:,2),R(:,2),par);
        Ts(:,1:2) = T_surface(:,1:2);

        Tf=zeros(par.nj,3);

        Tf = (0.0901 - 0.0575 * Sal) - 7.61e-4 * (par.g * h * par.rho_i / 1e4) + 273.16;

        del_h_i = h - h_init;
        del_h_o =  - (par.rho_i/par.rho_o) *  del_h_i;
        del_h_o_av = globmean(del_h_o(:,2),par);

        h_o_sal = h_o_s + del_h_o;

        h_o = h_o_init + globmean(del_h_o(:,2),par);

        ss = Sal_init;

        ss1 = Sal_init*par.h_o;

        Sal = (Sal_init*par.h_o)/h_o;

        v_sal = globmean(Sal(:,2),par).*h_o;

        iter_n = n;

        restartfolder_FIS = Directory+'/FISRestart';
        hhplot = hhplot(1:n+1);
        hhplotmax = hhplotmax(1:n+1);
        hhplotmin = hhplotmin(1:n+1);

        save(sprintf('%s/restart-exp-%.2d-Q-%.2d-1d-sphere-nonlinear-resnum-%.2d-eps-%.2d-CO2-%.2d-%s.mat',restartfolder_FIS,par.EBM_expnum,par.Qo,par.N,100*var1,par.df_CO2*10,var2),'h','R','Ta','To','Ts','T_f','Tf','Sal','h_o','qa','v_n','mask','theta','phi','iter_n','del_h_i','B','dhdt_cond','dhdt_odiff','hhplot','h_o','div_hv_np1','v_sal');

        if par.icelatlim==1 || par.icelatpole==1 || par.icelateq==1 || par.icethicklim==1 || par.icelatstable==1
            return;
        end

    end

    %% prepare for next time step:
    h(:,1)=h(:,2);
    h(:,2)=h(:,3);
    R(:,1)=R(:,2);
    R(:,2)=R(:,3);
    T_f(:,1)=T_f(:,2);
    T_f(:,2)=T_f(:,3);

    T_surface(:,1) = T_surface(:,2);
    T_ocean(:,1) = T_ocean(:,2);

    %% Reset h and R
    h(:,3) = 10000;
    R(:,3) = NaN;
    Sal(:,3)  = 50;

end
