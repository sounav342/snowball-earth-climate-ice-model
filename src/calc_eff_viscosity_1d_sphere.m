function [B,B_tN,B_tS,d_ftheta_v_N,d_ftheta_v_S]=calc_eff_viscosity_1d_sphere(~,par,T_f,T_surface,h,v,nn)
%% effective viscosity for 1d model in spherical coordinates
% nn is the exponent in Glen's Law
% T_surface is the surface air temp
% par.T_f is the freezing temp
% A_g is a function that calculates the temp dependence of ice viscosity, A(T)

% Comments
% There needs to be a line that says what happens if retreating ice ever equals par.Nedge = 1 or par.nj - XX

B=NaN(par.nj,1);
B_tN = NaN;
B_tS = NaN;
d_ftheta_v_N=NaN;
d_ftheta_v_S=NaN;
% BC on B - what if edge reaches pole?
B(1) = 0;
B(end) = 0;

dot_eps0=1e-14;
s = par.s(:);
c = par.c(:);
domain = par.domain(:);
interior = domain(domain ~= par.Nedge & domain ~= par.Sedge);

if ~isempty(interior)
    eps_xx = (v(interior) .* c(interior)) ./ (par.R .* s(interior));
    eps_yy = (v(interior + 1) - v(interior - 1)) ./ (2 * par.R * par.dtheta_rad);
    eps_zz = -(eps_xx + eps_yy);
    dot_eps = sqrt((eps_xx.^2 + eps_yy.^2 + eps_zz.^2) / 2);
    B(interior) = (h(interior) ./ par.R) ...
        .* exp(-0.03183 .* T_surface(interior) + 26.62) ...
        .* (dot_eps + dot_eps0).^(1 / nn - 1);
end

%% see ice_viscosity.m under AA folder in FISEBM/

% Solve for B at Nedge
if ~isnan(par.Nedge) && par.Nedge~=0 && any(domain==par.Nedge)
    j = par.Nedge;

    % Evaluate B_Tilde (@ j = par.Nedge)
    eps_xx = (1/(par.R*s(j)))*(v(j)*c(j));
    eps_yy = (1/par.R)*(v(j)-v(j-1))/(par.dtheta_rad);
    eps_zz = -(eps_xx + eps_yy);
    dot_eps=sqrt((eps_xx^2+eps_yy^2+eps_zz^2)/2);
    B_tN = (h(j)/par.R)*(exp(-0.03183*T_surface(j)+26.62))*(dot_eps+dot_eps0)^(1/nn-1);

    % Evaluate B_N (@ j = par.Nedge)
    d_ftheta_v_N = (1/4)*par.g*par.rho_i*h(j)^(2)*(1-par.mu)*B_tN^(-1) - (1/2)*v(j)*c(j)*s(j)^(-1);
    eps_xx = (1/(par.R*s(j)))*(v(j)*c(j));
    eps_yy = (1/par.R)*d_ftheta_v_N;
    eps_zz = -(eps_xx+eps_yy);
    dot_eps=sqrt((eps_xx^2+eps_yy^2+eps_zz^2)/2);

    B(j)=(h(j)/par.R)*(exp(-0.03183*T_surface(j)+26.62))*(dot_eps+dot_eps0)^(1/nn-1);
end

% Solve for B at Sedge
if ~isnan(par.Sedge) && par.Sedge~=0 && any(domain==par.Sedge)
    j = par.Sedge;

    % Evaluate B_tS (@ j = par.Sedge)
    eps_xx = (1/(par.R*s(j)))*(v(j)*c(j));
    eps_yy = (1/par.R)*(v(j+1)-v(j))/(par.dtheta_rad);
    eps_zz = -(eps_xx + eps_yy);
    dot_eps=sqrt((eps_xx^2+eps_yy^2+eps_zz^2)/2);

    B_tS = (h(j)/par.R)*(exp(-0.03183*T_surface(j)+26.62))*(dot_eps+dot_eps0)^(1/nn-1);

    % Evaluate B_tS (@ j = par.Sedge)
    d_ftheta_v_S = (1/4)*par.g*par.rho_i*h(j)^(2)*(1-par.mu)*B_tS^(-1) - (1/2)*v(j)*c(j)*s(j)^(-1);
    eps_xx = (1/(par.R*s(j)))*(v(j)*c(j));
    eps_yy = (1/par.R)*d_ftheta_v_S;
    eps_zz = -(eps_xx+eps_yy);
    dot_eps=sqrt((eps_xx^2+eps_yy^2+eps_zz^2)/2);
    B(j)=(h(j)/par.R)*(exp(-0.03183*T_surface(j)+26.62))*(dot_eps+dot_eps0)^(1/nn-1);
end
