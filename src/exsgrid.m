function [h,R] = exsgrid(h,par,R,S,exsgi,v_np1,n)
for i = 1:length(exsgi)
    j = exsgi(i);
    if isnan(R(j,3))
        % This code runs when sgi advance hasn't updated a R(eSC,3)
        if j <= par.EQ
            delta_h = par.dt*S(j);
            h_sg = R(j,2)*par.Hcr;

            if delta_h <= 0
                % Retreat of Ice Shelf
                if abs(delta_h) > h_sg
                    % subgrid cell j melts entirely
                    % additional melt categorized as heat input into the ocean
                    h(j,3) = 0;
                    R(j,3) = 0;
                    % apply remainder melt to ocean

                    % par.Nedge remains unchange
                elseif abs(delta_h) < h_sg
                    % subgrid cell j only partially melts
                    h(j,3) = 0;
                    h_sg = h_sg - abs(delta_h); % find new subgrid height
                    R(j,3) = find_R(h_sg,par); % update R in subgrid cell j
                end
            elseif delta_h > 0
                % Advance of Ice Shelf
                h_sg_2R1 = par.Hcr - h_sg;
                if delta_h > h_sg_2R1
                    % subgrid cells j fill entirely
                    h(j,3) = par.Hcr;
                    R(j,3) = 1;
                    % additional ice mass goes where?
                elseif delta_h < h_sg_2R1
                    % subgrid cell j partly fills
                    h(j,3) = 0;
                    h_sg = h_sg + delta_h; % partly fill subgrid
                    R(j,3) = find_R(h_sg,par); % update R at j
                end
            end
        elseif j > par.EQ
            delta_h = par.dt*S(j);
            h_sg = R(j,2)*par.Hcr;

            if delta_h <= 0
                % Retreat of Ice Shelf
                if abs(delta_h) > h_sg
                    % subgrid cell j melts entirely
                    % additional melt categorized as heat input into the ocean
                    h(j,3) = 0;
                    R(j,3) = 0;
                    % apply remainder melt to ocean:

                elseif abs(delta_h) < h_sg
                    % subgrid cell j only partially melts
                    h(j,3) = 0;
                    h_sg = h_sg - abs(delta_h); % find new subgrid height
                    R(j,3) = find_R(h_sg,par); % update R in subgrid cell j
                end
            elseif delta_h > 0
                % Advance of Ice Shelf
                h_sg_2R1 = par.Hcr - h_sg;
                if delta_h > h_sg_2R1
                    % subgrid cell j fills entirely
                    h(j,3) = par.Hcr;
                    R(j,3) = 1;
                    % additional ice mass goes where?
                elseif delta_h < h_sg_2R1
                    % subgrid cell j partially fills
                    h(j,3) = 0;
                    h_sg = h_sg + delta_h; % partially fill subgrid
                    R(j,3) = find_R(h_sg,par); % update R at j
                end
            end
        end
    elseif ~isnan(R(j,3))
        % This code runs when sgi advance has updated a R(eSC,3)
        % The sgrid function has already applied advection mass balance to
        % the eSC cell and updated it in R(j,3). Thus, h_sg uses
        % (and consequently updates) R(j,3), instead of R(j,2).
        % The source function is then applied to the cell.
        if j <= par.EQ
            delta_h = par.dt*S(j);

            if R(j,3)==1
                h_sg = h(j,3);
            else
                h_sg = R(j,3)*par.Hcr;
            end

            if delta_h <= 0
                % Retreat of Ice Shelf
                if abs(delta_h) > h_sg
                    % subgrid cell j melts entirely
                    % additional melt categorized as heat input into the ocean
                    h(j,3) = 0;
                    R(j,3) = 0;
                    % apply remainder melt to ocean

                    % par.Nedge remains unchange
                elseif abs(delta_h) < h_sg
                    % subgrid cell j only partially melts
                    h_sg = h_sg - abs(delta_h); % find new subgrid height
                    if h_sg >= par.Hcr
                        h(j,3) = h_sg;
                    else
                        h(j,3) = 0;
                    end
                    R(j,3) = find_R(h_sg,par); % update R in subgrid cell j
                end
            elseif delta_h > 0
                % Advance of Ice Shelf
                h(j,3) = h_sg + delta_h;
                if h(j,3)>=par.Hcr
                    % h stays the same
                    R(j,3) = find_R(h(j,3),par);
                else
                    % subgrid
                    R(j,3) = find_R(h(j,3),par);
                    h(j,3) = 0;
                end
            end
        elseif j > par.EQ
            delta_h = par.dt*S(j);

            if R(j,3)==1
                h_sg = h(j,3);
            else
                h_sg = R(j,3)*par.Hcr;
            end

            if delta_h <= 0
                % Retreat of Ice Shelf
                if abs(delta_h) > h_sg
                    % subgrid cell j melts entirely
                    % additional melt categorized as heat input into the ocean
                    h(j,3) = 0;
                    R(j,3) = 0;
                    % apply remainder melt to ocean:

                elseif abs(delta_h) < h_sg
                    % subgrid cell j only partially melts
                    h_sg = h_sg - abs(delta_h); % find new subgrid height
                    if h_sg >= par.Hcr
                        h(j,3) = h_sg;
                    else
                        h(j,3) = 0;
                    end
                    R(j,3) = find_R(h_sg,par); % update R in subgrid cell j
                end
            elseif delta_h > 0
                % Advance of Ice Shelf
                h(j,3) = h_sg + delta_h;
                if h(j,3)>=par.Hcr
                    % h stays the same
                    R(j,3) = find_R(h(j,3),par);
                else
                    % subgrid
                    R(j,3) = find_R(h(j,3),par);
                    h(j,3) = 0;
                end
            end
        end
    end
end
