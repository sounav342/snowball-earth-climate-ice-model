function FISEBM(var1,var2,var3,var4,var5,var6)
% Asynchronous Coupling of Floating Ice Sheets 1D Spherical Ice Flow Model
% with a version of the Pollard and Kasting, 2005 Energy Balance Model.
% This version of FISEBM is reading epsilon restarts and saving the plots
% in the FISEBM folder (/epsilon/FISEBM/).

set(0,'DefaultFigureVisible','off');

% internal ice flow clock - 10000 years

% Conditions for Aborting Model
icelatlim = 0;
icelatpole = 0;
icelateq = 0;
icelatstable = 0;
noice = 0;

% Time
year = 365*86400;
t_stable = 0;
t_stable_max = 1e4*year;
Nstart = var3;
Nend = 100;

CO2_forcing = var6;

expnum = var5;

Qo = var4;
Qo = Qo';
iceline_N = NaN(length(Qo),length(expnum));
iceline_S = NaN(length(Qo),length(expnum));

Directory = string(pwd);

% Create Folder to Export Plots
plotfolder = Directory+'/Figures';
if ~exist(sprintf('%s',plotfolder),'file')
    [~,~] = mkdir(sprintf('%s',plotfolder));
end

% Create Folder to Export Restart Files
restartfolder = Directory+'/FISEBMRestart';
if ~exist(sprintf('%s',restartfolder),'file')
    [~,~] = mkdir(sprintf('%s',restartfolder));
end

for k = 1:length(CO2_forcing)
    CO2_i = CO2_forcing(k);
    for i = 1:length(Qo)
        Q_i = Qo(i);
        for j = 1:length(expnum)
            EBM_expnum = expnum(j);

            for N = Nstart:Nend
                EBMpar = EBM(EBM_expnum,N,Q_i,var1,var2,CO2_i);

                FISpar = fis_1D_sphere(EBM_expnum,N,Q_i,var1,var2,CO2_i);

                if FISpar.noice==0 && FISpar.icelatpole==0

                    %%% Plotting FIS ouput over Time %%%
                    % Declare Initial Conditions
                    if N==1
                        Rplot_EBM = EBMpar.R_initial;
                        Rplot_FIS = FISpar.R_initial;

                        hplot_ave = EBMpar.h_init_ave;
                        hplot_max = max(EBMpar.h_initial);
                        hplot_min = min(EBMpar.h_initial);

                        Toplot_ave = EBMpar.To_init_ave;
                        Toplot_max = max(EBMpar.To_initial);
                        Toplot_min = min(EBMpar.To_initial);

                        Taplot_ave = FISpar.Ta_init_ave;
                        Taplot_max = max(FISpar.Ta_initial);
                        Taplot_min = min(FISpar.Ta_initial);

                        vplot_ave = 0;
                        vplot_max = 0;
                        vplot_min = 0;

                        nplot_FIS = 0;

                        Eplot_toa = 0;
                        Eplot_sfc = 0;

                        netEvap = 100;
                        netPrec = 100;
                        Aplot = 0.30;

                        tplot = 0;
                        masterclock = 0;
                    elseif N==Nstart && Nstart>1
                        r_filename = sprintf('%s/restart-FISEBM-expnum-%.2d-Q-%.2d-N-%.2d-eps-%.2d-CO2-%.2d-%s.mat',restartfolder,EBM_expnum,Q_i,N-1,100*var1,CO2_i,var2);
                        load(r_filename);
                        % all the above fields are inputed
                    end

                    dt = (FISpar.n/FISpar.nt)*FISpar.Time;
                    masterclock = masterclock + dt;

                    % Concatenate 1:N-1 to N
                    next_col = size(Rplot_EBM,2)+1;
                    Rplot_EBM(:,next_col) = EBMpar.R_end;
                    Rplot_FIS(:,next_col) = FISpar.R_end;

                    next_idx = numel(hplot_ave)+1;
                    hplot_ave(next_idx) = EBMpar.h_ave;
                    hplot_max(next_idx) = EBMpar.h_max;
                    hplot_min(next_idx) = EBMpar.h_min;

                    Toplot_ave(next_idx) = EBMpar.To_ave;
                    Toplot_max(next_idx) = EBMpar.To_max;
                    Toplot_min(next_idx) = EBMpar.To_min;

                    Taplot_ave(next_idx) = FISpar.Ta_ave;
                    Taplot_max(next_idx) = FISpar.Ta_max;
                    Taplot_min(next_idx) = FISpar.Ta_min;

                    vplot_ave(next_idx) = FISpar.v_ave;
                    vplot_max(next_idx) = FISpar.v_max;
                    vplot_min(next_idx) = FISpar.v_min;
                    nplot_FIS(next_idx) = FISpar.n;

                    Eplot_toa(next_idx) = EBMpar.Eplot_toa;
                    Eplot_sfc(next_idx) = EBMpar.Eplot_sfc;

                    netPrec(next_idx) = EBMpar.netPrec;
                    netEvap(next_idx) = EBMpar.netEvap;

                    Aplot(next_idx) = EBMpar.Aplot;

                    % Other Fields
                    tplot(next_idx) = masterclock;
                    x1 = tplot/(FISpar.year*1000); % time in 10^3 yrs

                    Nplot = [0:N];
                    x2 = Nplot;

                    m=5; n=1; pos = 2;
                    fig = figure('Visible','off');

                    % Number of time steps per ice flow run
                    subplot(m,n,pos-1);
                    plot(x2,nplot_FIS);
                    h3=title(sprintf('FIS #n per N'));

                    % Ice Thickness over Time
                    subplot(m,n,pos);
                    plot(x2,hplot_ave,x2,hplot_max,x2,hplot_min);
                    h3=title(sprintf('h'));
                    legend({'have','hmax','hmin'},'Position',[0.18 0.7 0.05 0.05],'FontSize',4);

                    % Ocean Temp over Time
                    subplot(m,n,pos+1);
                    plot(x1,Toplot_ave,x1,Toplot_max,x1,Toplot_min);
                    h3=title(sprintf('To'));
                    legend({'Toave','Tomax','Tomin'},'Position',[0.18 0.53 0.05 0.05],'FontSize',4);

                    % Mean Atmos Temp over Time
                    subplot(m,n,pos+2);
                    plot(x1,Taplot_ave,x1,Taplot_max,x1,Taplot_min);
                    h3=title(sprintf('Ta'));
                    legend({'Taave','Tamax','Tamin'},'Position',[0.18 0.355 0.05 0.05],'FontSize',4);

    %                 % Mean Atmos Temp over Time

                    subplot(m,n,pos+3);
                    plot(x1,vplot_ave);
                    h3=title('v (m/yr)');
                    h2 = xlabel('t (10^3)');
                    legend({'vave'},'Position',[0.18 0.18 0.05 0.05],'FontSize',4);

                    saveas(fig,sprintf('%s/hovmoller-Q-%.2d-exp-%.2d-eps-%.2d-CO2-%.2d-%s-part-1.png',plotfolder,Q_i,EBM_expnum,100*var1,CO2_i*10,var2));
                    close(fig);

                    m=2; n=1;
                    fig = figure('Visible','off');

                    % Hovmoller of R over time
                    subplot(m,n,1);
                    jb=2; je=FISpar.nj-1;
                    X = 90-FISpar.theta(je:-1:jb);
                    Y = x1;
                    Z = Rplot_EBM(je:-1:jb,:);
                    contourf(X,Y,Z','LineColor','none');
                    colormap(cool);
                    colorbar('FontSize',4,'Position',[0.92 0.83 0.02 0.1]);
                    h1=ylabel('t');
                    h3=title(sprintf('R - EBM'));

                    % Hovmoller of R over time
                    subplot(m,n,2);
                    jb=2; je=FISpar.nj-1;
                    X = 90-FISpar.theta(je:-1:jb);
                    Y = x1;
                    Z = Rplot_FIS(je:-1:jb,:);
                    contourf(X,Y,Z','LineColor','none');
                    colormap(cool);
                    colorbar('FontSize',4,'Position',[0.92 0.23 0.02 0.1]);
                    h1=ylabel('t');
                    h3=title(sprintf('R - FIS'));

                    saveas(fig,sprintf('%s/hovmoller-Q-%.2d-exp-%.2d-eps-%.2d-CO2-%.2d-%s-part-2.png',plotfolder,Q_i,EBM_expnum,100*var1,CO2_i*10,var2));
                    close(fig);

                    m=2; n=1;
                    fig = figure('Visible','off');

                    % ToA Energy Balance over N
                    subplot(m,n,1);
                    plot(x1,Eplot_toa);
                    h2=xlabel('t (10^3 yrs)');
                    h3=title(sprintf('ToA Energy Balance = %.2d',Eplot_toa(end)));

                    % Surface Energy Balance over N
                    subplot(m,n,2);
                    plot(x2,Eplot_sfc);
                    h2=xlabel('N');
                    h3=title(sprintf('Surface Energy Balance = %.2d',Eplot_sfc(end)));

                    saveas(fig,sprintf('%s/hovmoller-Q-%.2d-exp-%.2d-eps-%.2d-CO2-%.2d-%s-part-3.png',plotfolder,Q_i,EBM_expnum,100*var1,CO2_i*10,var2));
                    close(fig);

                    m=2; n=1;
                    fig = figure('Visible','off');

                    % Net Prec & Evap over N
                    subplot(m,n,1);
                    plot(x1,netPrec,x1,netEvap);
                    h2=xlabel('t (10^3 yrs)');
                    h3=title(sprintf('Net Prec = %.2d & Net Evap = %.2d',netPrec(end)*EBMpar.year*100/EBMpar.rho_i,netEvap(end)*EBMpar.year*100/EBMpar.rho_i));

                    % Global Mean Albedo
                    subplot(m,n,2);
                    plot(x2,Aplot);
                    h2=xlabel('N');
                    h3=title(sprintf('Global Mean Albedo = %.2d',Aplot(end)));

                    saveas(fig,sprintf('%s/hovmoller-Q-%.2d-exp-%.2d-eps-%.2d-CO2-%.2d-%s-part-4.png',plotfolder,Q_i,EBM_expnum,100*var1,CO2_i*10,var2));
                    close(fig);

                    % Iceline Stability reached after some time T
                end

                if FISpar.icelateq==1
                    FISpar.Nedge = FISpar.EQ;
                    FISpar.Sedge = FISpar.EQ;
                    iceline_N(i,j) = FISpar.Nedge;
                    iceline_S(i,j) = FISpar.Sedge;
                    fprintf(1,'Equatorial closure at N = %.2d: continuing coupled run as full snowball at N = %.2d.\n',N,N+1);
                end

                if FISpar.noice==1 || FISpar.icelatpole==1
                    if strcmp(FISpar.text,'start')==1
                        if FISpar.Nedge_initial==0
                            FISpar.Nedge_initial = FISpar.EQ;
                        end
                        if FISpar.Sedge_initial==0
                            FISpar.Sedge_initial = FISpar.EQ;
                        end
                        iceline_N(i,j) = FISpar.Nedge_initial;
                        iceline_S(i,j) = FISpar.Sedge_initial;
                        break;
                    elseif strcmp(FISpar.text,'end')==1
                        if FISpar.Nedge==0
                            FISpar.Nedge = FISpar.EQ;
                        end
                        if FISpar.Sedge==0
                            FISpar.Sedge = FISpar.EQ;
                        end
                        iceline_N(i,j) = FISpar.Nedge;
                        iceline_S(i,j) = FISpar.Sedge;
                        break;
                    end
                end

                save(sprintf('%s/restart-FISEBM-expnum-%.2d-Q-%.2d-N-%.2d-eps-%.2d-CO2-%.2d-%s.mat',restartfolder,EBM_expnum,Q_i,N,100*var1,CO2_i,var2),'iceline_N','iceline_S','Rplot_EBM','Rplot_FIS','hplot_ave','hplot_max','hplot_min','Toplot_ave','Toplot_max','Toplot_min','Taplot_ave','Taplot_max','Taplot_min','tplot','masterclock','vplot_ave','vplot_max','vplot_min','nplot_FIS','Eplot_toa','Eplot_sfc','netEvap','netPrec','Aplot');

                if N>=20 && FISpar.n==FISpar.nt
                    fprintf('Model Equilibrated: N>=55 & ice flow model reached 1Myr');
                    break;
                end

            end

            % Rewrite the Restart after each test
            if FISpar.noice==0 && FISpar.icelatpole==0
                save(sprintf('%s/restart-FISEBM-expnum-%.2d-Q-%.2d-N-%.2d-eps-%.2d-CO2-%.2d-%s.mat',restartfolder,EBM_expnum,Q_i,N,100*var1,CO2_i,var2),'iceline_N','iceline_S','Rplot_EBM','Rplot_FIS','hplot_ave','hplot_max','hplot_min','Toplot_ave','Toplot_max','Toplot_min','Taplot_ave','Taplot_max','Taplot_min','tplot','masterclock','vplot_ave','vplot_max','vplot_min','nplot_FIS','Eplot_toa','Eplot_sfc','netEvap','netPrec','Aplot');
            end

        end
    end
end

fprintf(1,'FISEBM, done.\n');
