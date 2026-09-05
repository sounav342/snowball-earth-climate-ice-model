function Pn = P(Ta,~,qa,par,~,~)
%P Calculate precipitation where relative humidity exceeds 0.8.

Pn = zeros(par.nj,1);

ra_max = 0.8;

rhoa = rho_a(Ta,par);
e = e_p(qa,par);
ew = e_w(Ta,par);
ew = mb_2_Pa(ew);
ra = r_a(e,ew);

index = ra > ra_max;

Pn(index) = (rhoa(index) .* par.h_q .* qa(index) ./ par.tau_p) ...
    .* (ra(index).^3 - ra_max^3);

if any(Pn(index)<0)
    Pn(Pn < 0) = 0;
    error('Snowball:NegativePrecipitation','Calculated precipitation is negative.');
end
end
