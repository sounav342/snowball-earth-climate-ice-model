function e=e_p(qa,par)
% partial pressure of water vapor (Pa)

x = par.epsilon^(-1)*qa./(1-qa);
e = par.p_a*x./(1+x);
