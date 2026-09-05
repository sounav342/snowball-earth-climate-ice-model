function [atm_E,ocn_E,netatm_E,netocn_E,atm_diff,ocn_diff,netatm_diff,netocn_diff] = hflux(par,Surf,ToA,diff_To_np1,diff_Ta_np1,diff_qa_np1)
% Calculates the total ocean and atmospheric heat transport
% as a function of latitude.
atm = zeros(par.nj,1);
ocn = zeros(par.nj,1);

% Heat Transport (based on energy imbalances)
s = par.s_col;
f = -(ToA(:) - Surf(:)).*s;
g = -Surf(:).*s;
phi = (1:(par.nj-1)).';
theta = par.theta(:);
theta_span = (theta(phi)-theta(1))*pi/180;
atm(1:par.nj-1) = theta_span .* cumsum(f(2:end)+f(1:end-1)) ./ (2*phi);
ocn(1:par.nj-1) = theta_span .* cumsum(g(2:end)+g(1:end-1)) ./ (2*phi);
mean_atm = specint(f,par.nj-1,par);
mean_ocn = specint(g,par.nj-1,par);

atm_E = (2*pi*par.R^2*(atm-mean_atm))/1e15; % output in PW
ocn_E= (2*pi*par.R^2*(ocn-mean_ocn))/1e15; % output in PW

netatm_E = specint(atm_E,par.EQ-1,par);
netocn_E = specint(ocn_E,par.EQ-1,par);

% Heat Transport (based on diffusion terms)

atm2 = zeros(par.nj,1);
ocn2 = zeros(par.nj,1);

% Heat Transport (based on diffusion terms)
f2 = diff_Ta_np1(:).*s + par.L_v*diff_qa_np1(:).*s;
g2 = diff_To_np1(:).*s;
atm2(1:par.nj-1) = theta_span .* cumsum(f2(2:end)+f2(1:end-1)) ./ (2*phi);
ocn2(1:par.nj-1) = theta_span .* cumsum(g2(2:end)+g2(1:end-1)) ./ (2*phi);
mean_atm2 = specint(f2,par.nj-1,par);
mean_ocn2 = specint(g2,par.nj-1,par);

atm_diff = 2*pi*par.R^2*(atm2-mean_atm2)/1e15;
ocn_diff = 2*pi*par.R^2*(ocn2-mean_ocn2)/1e15;

netatm_diff = 2*pi*par.R^2*(specint(atm_diff,par.EQ-1,par));
netocn_diff = 2*pi*par.R^2*(specint(ocn_diff,par.EQ-1,par));
