function Ts = Ts_extrap_FIS(h,R,Ts,par)
% This function extrapolates surface ice temperatures when ice
% cover advances to previously ice-free regions.
% It also returns the surface temperature field to NaN if
% ice-free regions are generated.

% Retreat to Open Ocean from Ice Cover
all_Ts_rt = R(:,2)==0 & R(:,1)>0;
Ts(all_Ts_rt) = NaN;

% General case of ice expanding to new latitudes
% i.e. advance from open ocean to ice cover

% Ts Extrapolation Function - 04/20/2020
newice = find(R(:,2)>0 & R(:,1)==0);

if ~isempty(newice)
    for i = 1:length(newice)
        j = newice(i);
        if j<par.EQ
            % Closest Ice GC Poleward of j in NH
            j_NH = find(R(1:j-1,2)>0 & R(1:j-1,1)>0,1,'last');
            Ts(j) = Ts(j_NH);
        elseif j>par.EQ
            % Closest Ice GC Poleward of j in SH
            j_SH = find(R((j+1):end,2)>0 & R((j+1):end,1)>0,1,'first') + j;
            Ts(j) = Ts(j_SH);
        else
            % Closest Ice GC Poleward of EQ
            j_NH = find(R(1:j-1,2)>0 & R(1:j-1,1)>0,1,'last');
            j_SH = find(R((j+1):end,2)>0 & R((j+1):end,1)>0,1,'first') + j;
            cNH = abs(j_NH-j);
            cSH = abs(j_SH-j);
            if min(cNH,cSH)==cNH
                % Ice cell in NH is closer
                Ts(j) = Ts(j_NH);
            elseif min(cNH,cSH)==cSH
                % Ice cell in SH is closer
                Ts(j) = Ts(j_SH);
            else
                % Same distance in either hemisphere
                Ts(j) = Ts(j_NH); % = Ts(j_SH)
            end
        end
    end
end

% Simplest Option: Simply set new ice to Ts=Tf
