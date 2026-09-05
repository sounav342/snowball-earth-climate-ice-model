function par = fis_1D_sphere(EBM_expnum,N,Qo,var1,var2,CO2_forcing)
%% Solve the model for Marine ice sheets
%% Eli, June 2011
%% Francisco, June 2017

par = set_FIS_parameters(EBM_expnum,N,Qo,CO2_forcing);

%% our 1d horizontal flow model using spherical coordinates:
par = integrate_h_1d_sphere(par,var1,var2);

fprintf(1,'done.\n');
