function [v,par]=solve_v_1d_matrix_form_sphere(n,par,B,B_tN,B_tS,h)
%% solve diagnostic momentum equation for 1d linearized model in
%% spherical coordinates, by writing it in matrix form.

%% initialize variables:
v=NaN(par.nj,1);
a=NaN(par.nj,1);
b=a; c=a; d=a;
s = par.s(:);
cc = par.c(:);

%% setup matrix and rhs:
domain = par.domain(:);
work_domain = domain(2:(end-1));
interior = work_domain(work_domain ~= par.Nedge & work_domain ~= par.Sedge);

if ~isempty(interior)
    B_jphalf = (B(interior) + B(interior + 1)) / 2;
    B_jmhalf = (B(interior) + B(interior - 1)) / 2;
    s_jphalf = (s(interior) + s(interior + 1)) / 2;
    s_jmhalf = (s(interior) + s(interior - 1)) / 2;

    a(interior) = B_jmhalf .* (s_jmhalf ./ s(interior) + s(interior - 1) ./ s_jmhalf) / par.dtheta_rad^2;
    b(interior) = -B_jphalf .* (s_jphalf ./ s(interior) + s(interior) ./ s_jphalf) / par.dtheta_rad^2 ...
        -B_jmhalf .* (s_jmhalf ./ s(interior) + s(interior) ./ s_jmhalf) / par.dtheta_rad^2 ...
        -(cc(interior) ./ s(interior)) .* B(interior) .* (cc(interior) ./ s(interior));
    c(interior) = B_jphalf .* (s_jphalf ./ s(interior) + s(interior + 1) ./ s_jphalf) / par.dtheta_rad^2;
    d(interior) = par.g * par.rho_i * (1 - par.mu) .* h(interior) .* (h(interior + 1) - h(interior - 1)) / (2 * par.dtheta_rad);
end

for i = 1:length(work_domain)
    j = work_domain(i);

    if j == par.Nedge
        % No calc of B_jphalf
        % No calc of s_jphalf
        % BC of c(j) = 0
        % New a(j), b(j), d(j)
        B_jphalf = NaN;
        B_jmhalf =(B(j)+B(j-1))/2;
        %% sine at half locations:
        s_jphalf = NaN;
        s_jmhalf =(s(j)+s(j-1))/2;
        %% v-equations:
        a(j) = s(j)^(-1)*B_jmhalf*s_jmhalf/par.dtheta_rad^2 + B_jmhalf*s_jmhalf^(-1)*s(j-1)/par.dtheta_rad^2;

        b(j) = B(j)*s(j)^(-1)*s(j+1)/par.dtheta_rad^2 - cc(j)*s(j)^(-1)*B(j)/(2*par.dtheta_rad) -cc(j)*s(j)^(-1)*B(j)*s(j)^(-1)*s(j+1)/(2*par.dtheta_rad) - B(j)/par.dtheta_rad^2 - B_jmhalf*s_jmhalf^(-1)*s(j)/par.dtheta_rad^2 - s(j)^(-1)*B_jmhalf*s_jmhalf/par.dtheta_rad^2 - cc(j)^(2)*s(j)^(-2)*B(j);

        c(j) = 0; % BC at terminus b/c v_Nedge+1 not defined
        test = 1; % check d(j) - artificially scaling back pressure - XX
        d(j) = par.g*par.rho_i*(1-par.mu)*h(j)*(h(j)-h(j-1))/(par.dtheta_rad) - par.g*par.rho_i*(1-par.mu*(test))*(h(j)^2)*B(j)*B_tN^(-1)*(1+s(j)^(-1)*s(j+1))/(4*par.dtheta_rad);

    elseif j == par.Sedge
        % No calc of B_jmhalf
        % No calc of s_jmhalf
        % BC of a(j) = 0
        % New a(j), b(j), d(j)
        B_jphalf = (B(j)+B(j+1))/2;
        B_jmhalf = NaN;
        %% sine at half locations:
        s_jphalf = (s(j)+s(j+1))/2;
        s_jmhalf = NaN;
        %% v-equations:
        a(j) = 0; % BC at terminus, v_Sedge-1 not defined

        b(j) = -s(j)^(-1)*B_jphalf*s_jphalf/par.dtheta_rad^2 - B_jphalf*s_jphalf^(-1)*s(j)/par.dtheta_rad^2 - B(j)/par.dtheta_rad^2 - s(j)^(-2)*cc(j)^(2)*B(j) + B(j)*s(j)^(-1)*s(j-1)/par.dtheta_rad^2 + B(j)*cc(j)*s(j)^(-1)/(2*par.dtheta_rad) + B(j)*cc(j)*s(j)^(-1)*s(j)^(-1)*s(j-1)/(2*par.dtheta_rad);

        c(j) = s(j)^(-1)*B_jphalf*s_jphalf/par.dtheta_rad^2 + B_jphalf*s_jphalf^(-1)*s(j+1)/par.dtheta_rad^2;
        test = 1; % check d(j) - artificially scaling back pressure - XX
        d(j) = par.g*par.rho_i*(1-par.mu)*h(j)*(h(j+1)-h(j))/(par.dtheta_rad) + par.g*par.rho_i*(1-par.mu*(test))*(h(j)^2)*B(j)*B_tS^(-1)*(1+s(j)^(-1)*s(j-1))/(4*par.dtheta_rad);
    end
end
test = 1; % XX
if par.Nedge~=0
    par.aN = a(par.Nedge);
    par.bN = b(par.Nedge);
    par.cN = c(par.Nedge); % c=0
    par.dN = d(par.Nedge);
    par.d1N = par.g*par.rho_i*(1-par.mu)*h(par.Nedge)*(h(par.Nedge)-h(par.Nedge-1))/(par.dtheta_rad);
    par.d2N = par.g*par.rho_i*(1-par.mu*(test))*(h(par.Nedge)^2)*B(par.Nedge)*B_tN^(-1)*(1+s(par.Nedge)^(-1)*s(par.Nedge+1))/(4*par.dtheta_rad);
end
if par.Sedge~=0
    par.aS = a(par.Sedge); % a=0
    par.bS = b(par.Sedge);
    par.cS = c(par.Sedge);
    par.dS = d(par.Sedge);
    par.d1S = par.g*par.rho_i*(1-par.mu)*h(par.Sedge)*(h(par.Sedge+1)-h(par.Sedge))/(par.dtheta_rad);
    par.d2S = par.g*par.rho_i*(1-par.mu*(test))*h(par.Sedge)^(2)*B(par.Sedge)*B_tS^(-1)*(1+s(par.Sedge)^(-1)*s(par.Sedge-1))/(4*par.dtheta_rad);
end

% Inspect a, b, c, d
%% prepare arrays for tridiagnonal solver:
%% note that v(1)=v(2)=v(par.nj-1)=v(par.nj-2)=0 and do not need to be solved for.

if ((par.Nedge==0) && (par.Sedge==0))
    % Solve Using Eli's Code
    A1=zeros(par.nj-4,1); B1=A1; C1=A1; D1=A1;
    A1(1)=0;
    B1(1)=b(3);
    C1(1)=c(3);
    D1(1)=d(3);
    A1(2:end-1)=a(4:par.nj-3);
    B1(2:end-1)=b(4:par.nj-3);
    C1(2:end-1)=c(4:par.nj-3);
    D1(2:end-1)=d(4:par.nj-3);
    A1(par.nj-4)=a(par.nj-2);
    B1(par.nj-4)=b(par.nj-2);
    C1(par.nj-4)=0;
    D1(par.nj-4)=d(par.nj-2);

    %% call tridiagonal solver:
    x = TDMAsolver(A1,B1,C1,D1);
    v(3:par.nj-2)=x;

    %% north and south b.c.:
    v(1)=0;
    v(2)=0;
    v(par.nj-1)=0;
    v(par.nj)=0;

elseif par.Sedge==0
    % Partial Ice NH, Total Ice SH
    A1N = zeros(par.Nedge-2,1); B1N=A1N; C1N=A1N; D1N=A1N;
    A1N(1)=0;
    B1N(1)=b(3);
    C1N(1)=c(3);
    D1N(1)=d(3);
    if par.Nedge >= 4
        A1N(2:end)=a(4:par.Nedge);
        B1N(2:end)=b(4:par.Nedge);
        C1N(2:end)=c(4:par.Nedge);
        D1N(2:end)=d(4:par.Nedge);
    end
    %% call tridiagonal solver for NH
    x = TDMAsolver(A1N,B1N,C1N,D1N);
    v(3:par.Nedge)=x;

    A1S = zeros(par.EQ-2,1); B1S=A1S; C1S=A1S; D1S=A1S;
    A1S(1:end-1)=a(par.EQ:par.nj-3);
    B1S(1:end-1)=b(par.EQ:par.nj-3);
    C1S(1:end-1)=c(par.EQ:par.nj-3);
    D1S(1:end-1)=d(par.EQ:par.nj-3);
    A1S(end)=a(par.nj-2);
    B1S(end)=b(par.nj-2);
    C1S(end)=0;
    D1S(end)=d(par.nj-2);
    %% call tridiagonal solver for SH
    x = TDMAsolver(A1S,B1S,C1S,D1S);
    v(par.EQ:par.nj-2)=x;

    % from array definition, v(par.Nedge+1:par.EQ-1) = NaN

    %% north pole b.c.:
    v(1)=0;
    v(2)=0;
    v(par.nj-1)=0;
    v(par.nj)=0;

elseif par.Nedge==0
    % Total Ice NH, Partial Ice SH
    A1N=zeros(par.EQ-2,1); B1N=A1N; C1N=A1N; D1N=A1N;
    A1N(1)=0;
    B1N(1)=b(3);
    C1N(1)=c(3);
    D1N(1)=d(3);
    A1N(2:end)=a(4:par.EQ);
    B1N(2:end)=b(4:par.EQ);
    C1N(2:end)=c(4:par.EQ);
    D1N(2:end)=d(4:par.EQ);
    %% call tridiagonal solver for the NH:
    x = TDMAsolver(A1N,B1N,C1N,D1N);
    v(3:par.EQ)=x;

    A1S=zeros(par.nj-par.Sedge-1,1); B1S=A1S; C1S=A1S; D1S=A1S;
    A1S(1:end-1)=a(par.Sedge:par.nj-3);
    B1S(1:end-1)=b(par.Sedge:par.nj-3);
    C1S(1:end-1)=c(par.Sedge:par.nj-3);
    D1S(1:end-1)=d(par.Sedge:par.nj-3);
    A1S(end)=a(par.nj-2);
    B1S(end)=b(par.nj-2);
    C1S(end)=0;
    D1S(end)=d(par.nj-2);

    %% call tridiagonal solver for the SH:
    x = TDMAsolver(A1S,B1S,C1S,D1S);
    v(par.Sedge:par.nj-2)=x;

    % from array definition, v(par.EQ+1:par.Sedge-1) = NaN
    %% north and south b.c.:
    v(1)=0;
    v(2)=0;
    v(par.nj-1)=0;
    v(par.nj)=0;

elseif ((par.Nedge~=0)&&(par.Sedge~=0))

    % Partial Ice NH, Partial Ice SH
    A1N = zeros(par.Nedge-2,1); B1N=A1N; C1N=A1N; D1N=A1N;
    A1N(1)=0;
    B1N(1)=b(3);
    C1N(1)=c(3);
    D1N(1)=d(3);
    if par.Nedge >= 4
        A1N(2:end)=a(4:par.Nedge);
        B1N(2:end)=b(4:par.Nedge);
        C1N(2:end)=c(4:par.Nedge);
        D1N(2:end)=d(4:par.Nedge);
    end
    %% call tridiagonal solver for NH
    x = TDMAsolver(A1N,B1N,C1N,D1N);
    v(3:par.Nedge)=x;

    A1S=zeros(par.nj-par.Sedge-1,1); B1S=A1S; C1S=A1S; D1S=A1S;
    A1S(1:end-1)=a(par.Sedge:par.nj-3);
    B1S(1:end-1)=b(par.Sedge:par.nj-3);
    C1S(1:end-1)=c(par.Sedge:par.nj-3);
    D1S(1:end-1)=d(par.Sedge:par.nj-3);
    A1S(end)=a(par.nj-2);
    B1S(end)=b(par.nj-2);
    C1S(end)=0;
    D1S(end)=d(par.nj-2);

    %% call tridiagonal solver for the SH:
    x = TDMAsolver(A1S,B1S,C1S,D1S);
    v(par.Sedge:par.nj-2)=x;

    %% north and south b.c.:
    v(1)=0;
    v(2)=0;
    v(par.nj-1)=0;
    v(par.nj)=0;
end
