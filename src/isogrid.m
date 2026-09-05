function [h,R] = isogrid(h,par,R,S,isogci,n)
% Mass Balance in Isolated Ice Shelves (grid cells interior to a subgrid cell)
% Grid to Subgrid Transitions carried out in main 'integrate' loop

for i = 1:length(isogci)
    j = isogci(i);

    delta_h = par.dt*S(j);
    if delta_h > 0
        % Ice Shelf Growth
        h(j,3) = h(j,2)+delta_h;
        R(j,3) = find_R(h(j,3),par);
    elseif delta_h <= 0
        % Ice Shelf Shrinkage
        h(j,3) = h(j,2) - abs(delta_h);

        if h(j,3)<0
            h(j,3) = 0;
            R(j,3) = 0;
        elseif h(j,3)<par.Hcr && h(j,3)>0
            R(j,3) = find_R(h(j,3),par);
            h(j,3) = 0;
        else
            R(j,3) = find_R(h(j,3),par);
            % h(j,3) as calculated above
        end

    end
end
