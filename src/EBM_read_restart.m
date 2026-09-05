function [h,R,Ta,To,Ts,T_f,qa,Sal,h_o,iter_n,del_h_i,h_o_sal,Sal_init]=EBM_read_restart(par,var2)
%% which_restart can be '1d-sphere-global','1d-sphere-partial','1d-sphere-icefree', '1d-sphere-modern'
h = zeros(par.nj,3);
To = h; Ta = h; qa = h; R=h; Sal=h;
Ts = NaN(par.nj,3);
T_f = ones(par.nj,3);
S = zeros(par.nj,1);
h_o = 4000;
iter_n=0;
del_h_i = zeros(par.nj,1);
h_o_sal = par.h_o*ones(par.nj,3);
Sal_init = 35*ones(par.nj,3);

which_restart = par.version;

Directory = string(pwd);

restartfolder_FIS = Directory+'/FISRestart';
restartfolder_EBM = Directory+'/EBMRestart';

if par.N==1 && par.EBM_expnum<=22
    fprintf(1,'Initializing EBM @ N = 1 with mkrestart file.\n');
    restart_filename = sprintf('%s/mkrestart-EBM-%s-Q-%.2d-eps-%.2d-CO2-%.2d-%s.mat',restartfolder_EBM,which_restart,par.Qo,100*par.A_a,par.df_CO2*10,var2);
else
    fprintf(1,'Initializing EBM @ N = %.2d with FIS restart file.\n',par.N);
    restart_filename = sprintf('%s/restart-exp-%.2d-Q-%.2d-1d-sphere-nonlinear-resnum-%.2d-eps-%.2d-CO2-%.2d-%s.mat',restartfolder_FIS,par.EBM_expnum,par.Qo,par.N-1,100*par.A_a,par.df_CO2*10,var2);
end

if exist(restart_filename,'file')
    load(restart_filename);
    Ts(:,3) = NaN;
    To(:,3) = 0;
else
    fprintf(2,'EBM error: no mkrestart or FIS restart file %s.\n',restart_filename);
    error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
end
