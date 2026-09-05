function ew=e_w(Ta,par)
% Saturation water vapor pressure (mb) from Emmanuel 4.4.14
% p 116-117
% over a planar surface of liquid water
TT = Ta - par.T_f;
ew = 6.112*exp((17.67*TT)./(TT+243.5));

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Saturation water vapor pressure (mb) from Emmanuel 4.4.14
% % p 116-117
% % over a planar surface of liquid water

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Saturation water vapor pressure (mb) from Emmanuel 4.4.14
% % p 116-117
% % over a planar surface of liquid water
%     % Ensure Ta and T_f are numerical arrays
%     % Calculate TT and ew
