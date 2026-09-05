function [h,R] = OOgrid(h,par,R,S,oogci)
% Mass Balance in Open Ocean Grid Cells
% Defined as an open ocean cell not directly adjacent to a shelf terminus
% Grid to Subgrid Transitions carried out in main 'integrate' loop
for i = 1:length(oogci)
    j = oogci(i);

    delta_h = par.dt*S(j);

    if delta_h > 0
        % Ice Growth on Ocean? XX
        if delta_h>par.Hcr
            h(j,3) = delta_h;
            R(j,3) = find_R(h(j,3),par);
        else
            R(j,3) = find_R(delta_h,par);
            h(j,3) = 0;
        end
    elseif delta_h <= 0
        % additional heat input to ocean? XX
        h(j,3) = 0;
        R(j,3) = 0;
    end
end
