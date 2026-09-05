function TS = Ta_surface_function(par)
% Calculate Initial Surface Air Temperature Profile
y_middle = par.y(par.nj-1)/2;
lat = (par.y-y_middle)*0.5*pi/y_middle;
if strcmp(par.Ta_surface_profile_type,'mine')
  %% Initial Try
  TS=par.T_f+(-20-40*sin(pi*(abs(par.y)-0.5*par.y(par.nj-1))/par.y(par.nj-1)).^4);
elseif strcmp(par.Ta_surface_profile_type,'cold')
  %% Dorian's low CO2:
  TS=par.T_f-79+48*cos(lat).^2;
elseif strcmp(par.Ta_surface_profile_type,'warm')
  %% Dorian's high CO2:
  TS=par.T_f-58+47*cos(lat).^2;
else
  disp('*** no such surface temperature profile.');
  quit
end

