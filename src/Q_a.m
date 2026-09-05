function Qa = Q_a(R, par)
% Qa: net solar radiative flux absorbed by atmosphere

% Define constants
Aa = 0.13;

Qa = Aa * par.Q/4 .* par.S_col;

end
