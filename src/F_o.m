function Fo=F_o(C,G,par,R)
% Net Freezing or Melting Rate at Ice Base in prognostic ice h eq
% C is conductive flux under ice, G is the geothermal heat flux onto ice
Fo = zeros(par.nj,1);
idx = R > 0;
Fo(idx) = (C(idx) - G(idx)) ./ (par.rho_i * par.L_f);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Net Freezing or Melting Rate at Ice Base in prognostic ice h eq
% % C is conductive flux under ice, G is the geothermal heat flux onto ice
% % Find indices where R is greater than 0
% % Calculate net freezing or melting rate using vectorized operations
