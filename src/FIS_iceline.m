function [par] = FIS_iceline(h,R,par,text)
% This function determines the ice line and ends FIS if:
% (0) no ice at initialization of FIS
% (1) ice retreats to the poles in FIS
% (2) ice first advances to the equator in FIS
% (3) ice line stabilizes
% (4) the global ice area changes more than the area
%     corresponding to an equatorial 4 deg lat strip
%     (enough to significantly feed back on climate)
% (5) global ice thickness changes enough to feed back on climate
%% Operation includes
% i) final round of plots
% ii) make restart file
% iii) FIS early termination

par.text = text;

%% Nedge/Sedge Update %%
if strcmp('start',text)==1
    t = 2;
elseif strcmp('end',text)==1
    t = 3;
end

arg_NH = find(h(1:par.EQ,t)==0,1,'first');
if ~isempty(arg_NH)
    if arg_NH==1
        Nedge=NaN;
    elseif arg_NH>=2
        Nedge=arg_NH-1;
    else
        fprintf(1,'error: edge\n');
        error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
    end
elseif sum(h(1:par.EQ,t)>0)==ceil(par.nj/2)
    % no h=0 found - ice covered
    Nedge=0;
else
    fprintf(1,'unforseen Nedge routine\n');
    error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
end

sh_len = length(h(par.EQ:end,t));
arg_SH = find(h(par.EQ:end,t)==0,1,'last');
if ~isempty(arg_SH)
    if arg_SH==sh_len
        Sedge=NaN;
    elseif arg_SH<=(par.nj-1)
        Sedge = arg_SH + par.EQ;
    else
        fprintf(1,'error: edge\n');
        error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
    end
elseif sum(h(par.EQ:end,t)>0)==ceil(par.nj/2)
    % no h=0 found - ice covered
    Sedge = 0;
else
    fprintf(1,'unforseen Nedge routine\n');
    error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
end

par.Nedge = Nedge;
par.Sedge = Sedge;

% Operation (0)
if (isnan(Nedge) || isnan(Sedge)) && t==3
    par.noice = 1;
    fprintf(1,sprintf('FISpar.noice = 1, EBM_expnum %.2d, n = %.2d\n',par.EBM_expnum,par.n));
elseif (isnan(Nedge) || isnan(Sedge)) && t==2
    par.noice = 1;
    fprintf(1,sprintf('FISpar.noice = 1, EBM_expnum %.2d, n = %.2d\n',par.EBM_expnum,0));
end

% Operation (1)
if (Nedge<=2 && Nedge>=1) || (Sedge>=(par.nj-1) && Sedge<=par.nj)
    if t==2
        par.icelatpole = 1;
        fprintf(1,sprintf('FISpar.icelatpole = 1, EBM_expnum %.2d, n = %.2d\n',par.EBM_expnum,1));
    elseif t==3
        par.icelatpole = 1;
        fprintf(1,sprintf('FISpar.icelatpole = 1, EBM_expnum %.2d, n = %.2d\n',par.EBM_expnum,par.n));
    end
end

% Operation (2): ice first reaches the equator.
% This ends the current partial-ice FIS segment and hands the newly global
% geometry back to EBM; the coupled model can then continue as a full snowball.
eq_has_ice = R(par.EQ,t)>0 || h(par.EQ,t)>0;
nh_at_eq = ~isnan(Nedge) && (Nedge==0 || Nedge>=par.EQ-1);
sh_at_eq = ~isnan(Sedge) && (Sedge==0 || Sedge<=par.EQ+1);
initial_eq_has_ice = par.R_initial(par.EQ)>0 || par.h_initial(par.EQ)>0;
if strcmp('end',text)==1 && ~initial_eq_has_ice && eq_has_ice && nh_at_eq && sh_at_eq
    par.icelateq = 1;
    fprintf(1,sprintf('FISpar.icelateq = 1, EBM_expnum %.2d, n = %.2d\n',par.EBM_expnum,par.n));
end

% Operation (3)
% final round of plots
% make restart file
% flag run for partial equilibrium
% update equilibrium plot

%%% Operation (4) %%%

if strcmp('end',text)==1 && par.icelateq==0

    % If global land ice area changes more than area of 4 degrees latitude strip
    xi = sum(par.R_initial);
    xf = sum(R(:,t));

    if abs(xi-xf)>=2
        par.icelatlim = 1;
        fprintf(1,sprintf('Operation 4 Early FIS Termination: FISpar.icelatlim = 1, EBM_expnum %.2d, FIS_n = %.2d\n',par.EBM_expnum,par.n));
        % flag run for early termination
    end
end

%%% Operation (5) %%%

if strcmp('end',text)==1 && par.icelateq==0
    full_snowball = all(R(:,t)>0);
    if full_snowball
        h_initial_ave = globmean(max(par.h_initial,0),par);
        dh_ave = globmean(abs(h(:,t)-par.h_initial),par);
        dh_limit = max(par.dh_ice_couple_abs,par.dh_ice_couple_frac*h_initial_ave);

    end
end
