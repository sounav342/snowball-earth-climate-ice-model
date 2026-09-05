function int = specint(f,phi,par)
% Take the integral from the NP up to a given colatitude
% This function is used in calculating the heat transport from E imbalances.
% phi is the index of the colatitude
f = f(:);
x = f(2:phi+1)+f(1:phi);
% check that par.theta is in latitude, not colatitude - xx
theta = par.theta(:);
int = ((theta(phi)-theta(1))*pi/180)*sum(x)/(2*phi);
