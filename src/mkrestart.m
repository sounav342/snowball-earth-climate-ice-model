function mkrestart(par,var2)
%MKRESTART Create the initial EBM restart for a configured experiment.

restartFolder = fullfile(pwd,'EBMRestart');
if ~exist(restartFolder,'dir')
    mkdir(restartFolder);
end

if par.EBM_expnum < 1 || par.EBM_expnum > 22 || par.EBM_expnum ~= fix(par.EBM_expnum)
    error('Snowball:InvalidExperiment','EBM_expnum must be an integer from 1 to 22.');
end

h = zeros(par.nj,3);
R = zeros(par.nj,3);
Ta = zeros(par.nj,3);
To = zeros(par.nj,3);
qa = zeros(par.nj,3);
Sal = ones(par.nj,3);
Ts = NaN(par.nj,3);
T_f = 273.16 .* ones(par.nj,3);

if par.EBM_expnum == 22
    [h,R,Ta,To,Ts,T_f,qa,Sal] = initialize_global_snowball( ...
        par,h,R,Ta,To,Ts,T_f,qa,Sal);
else
    [iceThickness,salinity,initializeOcean] = partial_parameters(par.EBM_expnum);
    [h,R,Ta,To,Ts,T_f,qa,Sal] = initialize_partial_cover( ...
        par,h,R,Ta,To,Ts,T_f,qa,Sal,iceThickness,salinity,initializeOcean);
end

restartFile = fullfile(restartFolder,sprintf( ...
    'mkrestart-EBM-%s-Q-%.2d-eps-%.2d-CO2-%.2d-%s.mat', ...
    char(par.version),par.Qo,100*par.A_a,par.df_CO2*10,var2));
save(restartFile,'h','qa','R','Ta','To','Ts','T_f','Sal');
end

function [iceThickness,salinity,initializeOcean] = partial_parameters(expnum)
iceThickness = 1000;
salinity = 35;
initializeOcean = true;

if expnum == 1
    iceThickness = 1;
elseif expnum == 3
    salinity = 50;
elseif expnum >= 4 && expnum <= 10
    iceThickness = 1500;
elseif expnum == 20
    initializeOcean = false;
end
end

function [h,R,Ta,To,Ts,T_f,qa,Sal] = initialize_partial_cover( ...
    par,h,R,Ta,To,Ts,T_f,qa,Sal,iceThickness,salinity,initializeOcean)

fullIce = [1:par.Nedge,par.Sedge:par.nj];
partialIce = [(par.Nedge+1):(par.Nedge+10),(par.Sedge-10):(par.Sedge-1)];
surfaceIce = [1:(par.Nedge+10),(par.Sedge-10):par.nj];

h(fullIce,:) = iceThickness;
Sal(:,:) = salinity;
if initializeOcean
    T_f = (0.0901 - 0.0575 .* Sal) ...
        - 7.61e-4 .* (par.g .* h .* par.rho_i ./ 1e4) + 273.16;
    To = T_f;
end

Ta(:,3) = Ta_surface_function(par);
Ta(:,1:2) = Ta(:,[3 3]);
Ts(surfaceIce,1:2) = Ta(surfaceIce,1:2) - 10;
R(fullIce,:) = 1;
R(partialIce,:) = 0.5;
end

function [h,R,Ta,To,Ts,T_f,qa,Sal] = initialize_global_snowball( ...
    par,h,R,Ta,To,Ts,T_f,qa,Sal)

h(:,:) = 1500;
Sal(:,:) = 50;
T_f = (0.0901 - 0.0575 .* Sal) ...
    - 7.61e-4 .* (par.g .* h .* par.rho_i ./ 1e4) + 273.16;
To = T_f;
Ta(:,3) = Ta_surface_function(par);
Ta(:,1:2) = Ta(:,[3 3]);
Ts = Ta;
Ta(:,:) = 210;
R(:,:) = 1;
end
