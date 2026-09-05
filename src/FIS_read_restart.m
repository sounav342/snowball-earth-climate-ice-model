function [h,R,Ta,To,Ts,T_f,qa,Hcr,S,Sal,h_o,h_o_sal,Sal_init]=FIS_read_restart(par,var1,var2)

B=NaN;

Directory = string(pwd);

% Restart Filename Associated with previous EBM Run
restartfolder_EBM = Directory+'/EBMRestart';

restart_filename = sprintf('%s/restart-EBM-expnum-%.2d-Q-%.2d-resnum-%.2d-eps-%.2d-CO2-%.2d-%s.mat',restartfolder_EBM,par.EBM_expnum,par.Qo,par.N,100*var1,par.df_CO2*10,var2);

if exist(restart_filename,'file')
    fprintf(1,'initializing FIS @ N = %.2d with restart from the last EBM run.\n',par.N);
    B=ones(1,par.nj)*1e16*par.h0/par.R; % xx
    load(restart_filename);

    % Set new fields to unusual defaults to catch bugs
    h(:,3) = NaN;
    R(:,3) = NaN;

else
    fprintf(1,'FIS error: no EBM restart file %s.\n',restart_filename);
    error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
end
