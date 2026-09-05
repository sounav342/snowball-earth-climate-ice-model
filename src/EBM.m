function par = EBM(EBM_expnum,N,Qo,var1,var2,CO2_forcing)
%% Solve the energy balance model for ocean, atmosphere, and ice
%% Francisco, 2018

par=set_EBM_parameters(EBM_expnum,N,Qo,var1,var2,CO2_forcing);

if strcmp(par.model_to_run,'2d-sphere')
    %% our 2D energy balance model using spherical coordinates:
    par = integrate_EBM_1d_sphere(par,var2);
else
    disp('*** no such model to run.');
    return
end

fprintf(1,'done.\n');
