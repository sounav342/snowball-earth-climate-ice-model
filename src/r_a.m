function ra=r_a(e,ew)
% atmospheric relative humidity
% all inputs must be in Pa
% ew_o output in mb, so convert to Pa first
% ew_s output in Pa
ra = (e./ew);
